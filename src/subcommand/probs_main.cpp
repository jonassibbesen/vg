/** \file probs_main.cpp
 *
 * Defines the "vg probs" subcommand.
 */

#include <unistd.h>
#include <getopt.h>
#include <chrono>

#include "subcommand.hpp"

#include "../gbwt_helper.hpp"
#include "../types.hpp"
 
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const double indel_prob = prob_to_logprob(0.0001);
const bool use_score = false; 

const double frag_length_mean = 277;
const double frag_length_sd = 43;

typedef pair<vector<gbwt::size_type>, int32_t> fragment_path_t;

//#define debug


pair<unordered_map<uint32_t, uint32_t>, vector<vector<uint32_t> > > find_connected_components(const unordered_map<int32_t, unordered_set<int32_t> > & connected_paths, const int32_t num_paths) {

    unordered_map<uint32_t, uint32_t> connected_path_components_1;
    vector<vector<uint32_t> > connected_path_components_2;

    for (uint32_t i = 0; i < num_paths; ++i) {

        if (connected_path_components_1.count(i) == 0) {

            std::queue<int32_t> search_queue;
            search_queue.push(i);

            connected_path_components_2.emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto connected_path_components_1_it = connected_path_components_1.emplace(search_queue.front(), connected_path_components_2.size() - 1);
                assert(connected_path_components_1_it.first->second == connected_path_components_2.size() - 1);
                
                if (connected_path_components_1_it.second) {

                    connected_path_components_2.back().emplace_back(search_queue.front());
                    auto out_edges_it = connected_paths.find(search_queue.front());

                    if (out_edges_it != connected_paths.end()) {

                        for (auto & out_edge: out_edges_it->second) {

                            if (connected_path_components_1.count(out_edge) == 0) {

                                search_queue.push(out_edge);
                            }
                        }
                    }
                }

                search_queue.pop();
            }

            sort(connected_path_components_2.back().begin(), connected_path_components_2.back().end());
        }
    }

    return make_pair(connected_path_components_1, connected_path_components_2);
}

double calc_read_prob(const Alignment & align, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs) {

    double read_path_prob = 0;

    auto & base_qualities = align.quality();
    int32_t cur_pos = 0;

    for (auto & mapping: align.path().mapping()) {

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

pair<gbwt::SearchState, int32_t> get_read_paths_search(const Path & read_path, const gbwt::GBWT & paths_index) {
            
    auto mapping_it = read_path.mapping().cbegin();
    
    gbwt::SearchState read_path_search = paths_index.find(mapping_to_gbwt(*mapping_it));
    int32_t read_path_length = mapping_to_length(*mapping_it);
    
    ++mapping_it;

    while (mapping_it != read_path.mapping().cend()) {

        read_path_search = paths_index.extend(read_path_search, mapping_to_gbwt(*mapping_it));
        read_path_length += mapping_to_length(*mapping_it);
        ++mapping_it;
    }

    return make_pair(read_path_search, read_path_length);
}

pair<vector<gbwt::size_type>, int32_t> get_read_path_ids(const Path & read_path, const gbwt::GBWT & paths_index) {

    auto read_paths_search = get_read_paths_search(read_path, paths_index);
    return make_pair(paths_index.locate(read_paths_search.first), read_paths_search.second);
}

void find_fragment_paths(vector<pair<vector<gbwt::size_type>, int32_t> > * fragment_paths, const pair<gbwt::SearchState, int32_t> & start_search, const Path & end_read_path, const gbwt::GBWT & paths_index, const xg::XG & xg_index, const int32_t max_pair_distance) {

    assert(!start_search.first.empty());

    std::queue<pair<gbwt::SearchState, int32_t> > fragment_paths_queue;
    fragment_paths_queue.push(start_search);

    assert(end_read_path.mapping_size() > 0);
    auto end_read_path_start = mapping_to_gbwt(end_read_path.mapping(0));

    // Perform depth-first path extension.
    while (!fragment_paths_queue.empty()) {

        auto & cur_fragment_path_search = fragment_paths_queue.front();

        // 
        if (cur_fragment_path_search.second > max_pair_distance) {

            fragment_paths_queue.pop();
            continue;                
        }

        // Stop current extension if end node is reached.
        if (cur_fragment_path_search.first.node == end_read_path_start) {

            auto mapping_it = end_read_path.mapping().cbegin();
            ++mapping_it;

            while (mapping_it != end_read_path.mapping().cend()) {

                cur_fragment_path_search.second += xg_index.node_length(mapping_it->position().node_id()); 
                cur_fragment_path_search.first = paths_index.extend(cur_fragment_path_search.first, mapping_to_gbwt(*mapping_it));
                ++mapping_it;
            }

            if (!cur_fragment_path_search.first.empty() and cur_fragment_path_search.second <= max_pair_distance) {
    
                fragment_paths->emplace_back(paths_index.locate(cur_fragment_path_search.first), cur_fragment_path_search.second);
            }

            fragment_paths_queue.pop();
            continue;            
        }
        
        auto out_edges = paths_index.edges(cur_fragment_path_search.first.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            fragment_paths_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 
        ++out_edges_it;

        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_search = paths_index.extend(cur_fragment_path_search.first, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_search.empty()) { 

                    fragment_paths_queue.push(make_pair(extended_search, cur_fragment_path_search.second + xg_index.node_length(gbwt::Node::id(out_edges_it->first))));
                }
            }

            ++out_edges_it;
        }

        if (out_edges.begin()->first != gbwt::ENDMARKER) {

            cur_fragment_path_search.second += xg_index.node_length(gbwt::Node::id(out_edges.begin()->first)); 
            cur_fragment_path_search.first = paths_index.extend(cur_fragment_path_search.first, out_edges.begin()->first);
        
            // End current extension if empty (no path found). 
            if (cur_fragment_path_search.first.empty()) { fragment_paths_queue.pop(); }
        
        } else {

            fragment_paths_queue.pop();
        }
    }
}

vector<fragment_path_t> get_fragment_paths_ids(const Path & read_path_1, const Path & read_path_2, const gbwt::GBWT & paths_index, const xg::XG & xg_index, const int32_t max_pair_distance) {

    vector<pair<vector<gbwt::size_type>, int32_t> > fragment_paths;

    function<size_t(const int64_t)> node_length_func = [&xg_index](const int64_t node_id) { return xg_index.node_length(node_id); };

    auto read_paths_1 = get_read_paths_search(read_path_1, paths_index);
    if (!read_paths_1.first.empty()) {

        const Path read_path_2_rc = reverse_complement_path(read_path_2, node_length_func);
        find_fragment_paths(&fragment_paths, read_paths_1, read_path_2_rc, paths_index, xg_index, max_pair_distance);
    }

    auto read_paths_2 = get_read_paths_search(read_path_2, paths_index);
    if (!read_paths_2.first.empty()) {

        const Path read_path_1_rc = reverse_complement_path(read_path_1, node_length_func);
        find_fragment_paths(&fragment_paths, read_paths_2, read_path_1_rc, paths_index, xg_index, max_pair_distance);
    }

    return fragment_paths;
}

void help_probs(char** argv) {
    cerr << "\nusage: " << argv[0] << " probs [options] <graph.xg> <paths.gbwt> <reads.gam(p)> > read_path_probs.txt" << endl
         << "options:" << endl
         << "    -t, --threads INT          number of compute threads to use" << endl
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
                {"threads",  no_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "t:ph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 't':
            num_threads = parse<int>(optarg);
            omp_set_num_threads(num_threads);
            break;
            
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

    assert(num_threads > 0);

    double time1 = gcsa::readTimer();

    unique_ptr<xg::XG> xg_index = vg::io::VPKG::load_one<xg::XG>(get_input_file_name(optind, argc, argv));

    double time2 = gcsa::readTimer();
    cerr << "Load XG " << time2 - time1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    unique_ptr<gbwt::GBWT> paths_index = vg::io::VPKG::load_one<gbwt::GBWT>(get_input_file_name(optind, argc, argv));

    double time3 = gcsa::readTimer();
    cerr << "Load GBWT " << time3 - time2 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

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

    function<size_t(const int64_t)> node_length_func = [&xg_index](const int64_t node_id) { return xg_index->node_length(node_id); };

    int num_fragments = 0;

    vector<unordered_map<int32_t, unordered_set<int32_t> > > connected_paths_threads(num_threads);
    vector<vector<vector<fragment_path_t> > > all_fragment_paths_threads(num_threads);

    double time4 = gcsa::readTimer();
    cerr << "Init structures " << time4 - time3 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    get_input_file(optind, argc, argv, [&](istream& in) {

        vg::io::for_each_interleaved_pair_parallel<Alignment>(in, [&](Alignment& align_1, Alignment& align_2) {

            if (align_1.has_path() and align_2.has_path()) {

#ifdef debug
                cerr << pb2json(align_1) << "\n";
                cerr << pb2json(align_2) << "\n";        

                auto read_paths_1_fw = get_read_path_ids(align_1.path(), *paths_index);
                auto read_paths_1_rc = get_read_path_ids(reverse_complement_path(align_1.path(), node_length_func), *paths_index);
                auto read_paths_2_fw = get_read_path_ids(align_2.path(), *paths_index);
                auto read_paths_2_rc = get_read_path_ids(reverse_complement_path(align_2.path(), node_length_func), *paths_index);

                for (auto & bla: read_paths_1_fw.first) {
                    cerr << bla << " ";
                }
                cerr << "| " << read_paths_1_fw.second << endl; 
                for (auto & bla: read_paths_1_rc.first) {
                    cerr << bla << " ";
                }
                cerr << "| " << read_paths_1_rc.second << endl; 

                for (auto & bla: read_paths_2_fw.first) {
                    cerr << bla << " ";
                }
                cerr << "| " << read_paths_2_fw.second << endl; 
                for (auto & bla: read_paths_2_rc.first) {
                    cerr << bla << " ";
                }
                cerr << "| " << read_paths_2_rc.second << endl; 
#endif 

                auto fragment_paths = get_fragment_paths_ids(align_1.path(), align_2.path(), *paths_index, *xg_index, frag_length_mean + 10 * frag_length_sd);

#ifdef debug            
                for (auto & bla: fragment_paths) {
                    for (auto & bla2: bla.first) {
                        cerr << bla2 << " ";
                    }
                    cerr << "| " << bla.second << " | " << normal_pdf<float>(bla.second, frag_length_mean, frag_length_sd) << endl; 
                }
#endif 

                if (!fragment_paths.empty()) {

                    auto anchor_path = fragment_paths.front().first.front();

                    for (auto & fragment_path: fragment_paths) {

                        for (uint32_t i = 0; i < fragment_path.first.size(); ++i) {

                            if (anchor_path != fragment_path.first.at(i)) {

                                connected_paths_threads.at(omp_get_thread_num())[anchor_path].emplace(fragment_path.first.at(i));
                                connected_paths_threads.at(omp_get_thread_num())[fragment_path.first.at(i)].emplace(anchor_path);
                            }
                        }
                    }

                    all_fragment_paths_threads.at(omp_get_thread_num()).emplace_back(fragment_paths);

                    if (all_fragment_paths_threads.at(omp_get_thread_num()).size() % 100000 == 0) {

                        cerr << omp_get_thread_num() << ": " << all_fragment_paths_threads.at(omp_get_thread_num()).size() << endl;        
                    }
                }

                // auto read_prob_1 = calc_read_prob(align_1, quality_match_probs, quality_mismatch_probs);
                // auto read_prob_2 = calc_read_prob(align_2, quality_match_probs, quality_mismatch_probs);
            }
        });
    });

    double time5 = gcsa::readTimer();
    cerr << "Found paths " << time5 - time4 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    for (uint32_t i = 1; i < connected_paths_threads.size(); ++i) {

        for (auto & connected_paths: connected_paths_threads.at(i)) {

            auto connected_paths_threads_it = connected_paths_threads.front().emplace(connected_paths.first, unordered_set<int32_t>());
            for (auto & paths: connected_paths.second) {

                connected_paths_threads_it.first->second.emplace(paths);
            }
        }
    }

    double time6 = gcsa::readTimer();
    cerr << "Merged threads " << time6 - time5 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    auto connected_path_components = find_connected_components(connected_paths_threads.front(), num_paths);

    cerr << connected_path_components.first.size() << endl;
    cerr << connected_path_components.second.size() << endl;

    double time62 = gcsa::readTimer();
    cerr << "Found components " << time62 - time6 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
 
    vector<vector<vector<fragment_path_t> > > clustered_fragment_paths(connected_path_components.second.size());

    for (auto & all_fragment_paths: all_fragment_paths_threads) {

        for (auto & fragment_paths: all_fragment_paths) {

            clustered_fragment_paths.at(connected_path_components.first.at(fragment_paths.front().first.front())).push_back(move(fragment_paths));
        }
    }

    double time7 = gcsa::readTimer();
    cerr << "Cluster fragment " << time7 - time6 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    for (uint32_t i = 0; i < clustered_fragment_paths.size(); ++i) {

        cout << "#";
        for (auto & connected_path: connected_path_components.second.at(i)) {

            cout << " " << connected_path;
        }
        cout << endl;
    }
 
    double time8 = gcsa::readTimer();
    cerr << "Wrote out " << time8 - time7 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    return 0;
}

// Register subcommand
static Subcommand vg_probs("probs", "calculate read path probabilities", main_probs);

