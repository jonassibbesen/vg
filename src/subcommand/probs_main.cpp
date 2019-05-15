/** \file probs_main.cpp
 *
 * Defines the "vg probs" subcommand.
 */

#include <unistd.h>
#include <getopt.h>
#include <chrono>

#include "subcommand.hpp"

#include "../gbwt_helper.hpp"
 
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const double indel_prob = prob_to_logprob(0.0001);


double calc_read_prob(const Alignment & alignment, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs) {

    double read_path_prob = 0;

    auto & base_qualities = alignment.quality();
    int32_t cur_pos = 0;

    for (auto & mapping: alignment.path().mapping()) {

        for (auto & edit: mapping.edit()) {

            if (edit_is_match(edit)) {

                for (int32_t i = cur_pos; i < cur_pos + edit.from_length(); ++i) {

                    read_path_prob += quality_match_probs.at(int32_t(base_qualities.at(i)));
                }

            } else if (edit_is_sub(edit)) {

                for (int32_t i = cur_pos; i < cur_pos + edit.from_length(); ++i) {

                    read_path_prob += quality_mismatch_probs.at(int32_t(base_qualities.at(i)));
                }
            
            } else if (edit_is_insertion(edit)) {

                read_path_prob += edit.to_length() * indel_prob;

            } else if (edit_is_deletion(edit)) {

                read_path_prob += edit.from_length() * indel_prob;
            }
        }
    } 

    return read_path_prob;
}

void help_probs(char** argv) {
    cerr << "\nusage: " << argv[0] << " probs [options] <graph.xg> <paths.gbwt> <reads.gam(p)> > read_path_probs.txt" << endl
         << "options:" << endl
         << "    -p, --progress             show progress" << endl
         << "    -h, --help                 print help message" << endl
         << endl;
}

int32_t main_probs(int32_t argc, char** argv) {

    if (argc == 2) {
        help_probs(argv);
        return 1;
    }
    
    int32_t num_threads = 1;
    bool show_progress = false;

    int32_t c;
    optind = 2;

    while (true) {
        static struct option long_options[] =
            {
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "ph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {


        case 'p':
            show_progress = true;
            break;

        case 'h':
        case '?':
            help_probs(argv);
            exit(1);
            break;

        default:
            abort();
        }
    }

    if (argc < optind + 1) {
        help_probs(argv);
        return 1;
    }

    double first = gcsa::readTimer();

    unique_ptr<xg::XG> xg_index = vg::io::VPKG::load_one<xg::XG>(get_input_file_name(optind, argc, argv));

    double second = gcsa::readTimer();
    cerr << "Load XG " << second - first << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    unique_ptr<gbwt::GBWT> paths_index = vg::io::VPKG::load_one<gbwt::GBWT>(get_input_file_name(optind, argc, argv));

    double third = gcsa::readTimer();
    cerr << "Load GBWT " << third - second << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    auto num_paths = paths_index->metadata.haplotype_count;
    cerr << num_paths << endl;

    vector<int32_t> out_data;

    vector<double> quality_match_probs;
    quality_match_probs.reserve(61);

    vector<double> quality_mismatch_probs;
    quality_mismatch_probs.reserve(61);

    for (double i = 0; i < 61; ++i) {

        quality_match_probs.emplace_back(prob_to_logprob(1 - pow(10, -1*i/10)));
        quality_mismatch_probs.emplace_back(prob_to_logprob(pow(10, -1*i/10) / 3));
    }

    function<void(Alignment &)> lambda = [&](Alignment & alignment) {

        cerr << pb2json(alignment) << "\n";

        const Path & path = alignment.path();

        gbwt::vector_type path_nodes;
        path_nodes.reserve(path.mapping_size());

        for (auto & mapping: path.mapping()) {

            path_nodes.emplace_back(mapping_to_gbwt(mapping));
        }

        gbwt::SearchState search_state = paths_index->find(path_nodes.begin(), path_nodes.end());
        auto path_names_fwd = paths_index->locate(search_state);

        auto path_rv = reverse_complement_path(path, [&xg_index](const int64_t node_id) { return xg_index->node_length(node_id); });
        path_nodes.clear();

        for (auto & mapping: path_rv.mapping()) {

            path_nodes.emplace_back(mapping_to_gbwt(mapping));
        }

        search_state = paths_index->find(path_nodes.begin(), path_nodes.end());
        auto path_names_rv = paths_index->locate(search_state);

        cerr << path_names_fwd.size() << " " << path_names_rv.size() << endl;

        path_names_fwd.insert(path_names_fwd.end(), path_names_rv.begin(), path_names_rv.end());

        for (auto & bla: path_names_fwd) {

            cerr << bla << " ";
        }

        cerr << endl;  

        cerr << calc_read_prob(alignment, quality_match_probs, quality_mismatch_probs) << endl;
    };


    get_input_file(optind, argc, argv, [&](istream& in) {
        vg::io::for_each(in, lambda);
    });

    return 0;
}

// Register subcommand
static Subcommand vg_probs("probs", "calculate read path probabilities", main_probs);

