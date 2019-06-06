/** \file probs_main.cpp
 *
 * Defines the "vg probs" subcommand.
 */

#include <unistd.h>
#include <getopt.h>
#include <chrono>
#include <limits>

#include "subcommand.hpp"

#include "../gbwt_helper.hpp"
#include "../types.hpp"
#include "../alignment.hpp"
 
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

//#define debug

const double frag_length_mean = 277;
const double frag_length_sd = 43;

static const double double_precision = numeric_limits<double>::epsilon() * 100;


struct AlignmentPath {

    gbwt::SearchState path; 
    vector<gbwt::size_type> path_ids;
    
    int32_t length;

    pair<int32_t, int32_t> scores;
    pair<int32_t, int32_t> mapqs;

    AlignmentPath() : length(0) {

        scores = make_pair(0,0);
        mapqs = make_pair(0,0);
    }
};

bool unordered_pair_compare(const pair<int32_t, int32_t> & pair_1, const pair<int32_t, int32_t> & pair_2) {

    if (pair_1.first == pair_2.first && pair_1.second == pair_2.second) {

        return true;
    
    } else if (pair_1.first == pair_2.second && pair_1.second == pair_2.first) {

        return true;
    } 

    return false;
}

bool operator==(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return (lhs.path_ids == rhs.path_ids && lhs.length == rhs.length && unordered_pair_compare(lhs.scores, rhs.scores) && unordered_pair_compare(lhs.mapqs, rhs.mapqs));
}

bool operator!=(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    return !(lhs == rhs);
}

bool operator<(const AlignmentPath & lhs, const AlignmentPath & rhs) { 

    if (lhs.path_ids.size() != rhs.path_ids.size()) {

        return (lhs.path_ids.size() < rhs.path_ids.size());    
    } 

    for (size_t i = 0; i < lhs.path_ids.size(); ++i) {

        if (lhs.path_ids.at(i) != rhs.path_ids.at(i)) {

            return (lhs.path_ids.at(i) < rhs.path_ids.at(i));    
        }         
    }

    if (!unordered_pair_compare(lhs.mapqs, rhs.mapqs)) {

        return (lhs.mapqs < rhs.mapqs);
    }
 
    if (lhs.length != rhs.length) {

        return (lhs.length < rhs.length);
    }

    if (!unordered_pair_compare(lhs.scores, rhs.scores)) {

        return (lhs.scores < rhs.scores);
    }

    return false;
}

bool AlignmentPathsSorter(const vector<AlignmentPath> & align_paths_1, const vector<AlignmentPath> & align_paths_2) {

    if (align_paths_1.size() != align_paths_2.size()) {

        return (align_paths_1.size() < align_paths_2.size());    
    } 

    if (align_paths_1.size() == 1) {

        return (align_paths_1.front() < align_paths_2.front());            
    }

    for (size_t i = 0; i < align_paths_1.size(); ++i) {

        if (align_paths_1.at(i) != align_paths_2.at(i)) {

            return (align_paths_1.at(i) < align_paths_2.at(i));    
        }         
    }   

    return false;     
}

ostream& operator<<(ostream& os, const AlignmentPath & align_path) {

    for (auto & id: align_path.path_ids) {

        os << id << " ";
    }

    os << "| " << align_path.length;
    os << " | (" << align_path.scores.first << ", " << align_path.scores.second << ")";
    os << " | (" << align_path.mapqs.first << ", " << align_path.mapqs.second << ")";

    return os;
}

struct PathClusters {

    vector<uint32_t> path_to_cluster_index;
    vector<vector<uint32_t> > cluster_to_path_index;
};

PathClusters find_path_clusters(const unordered_map<int32_t, unordered_set<int32_t> > & connected_paths, const int32_t num_paths) {

    PathClusters path_clusters;
    path_clusters.path_to_cluster_index = vector<uint32_t>(num_paths, -1);

    for (size_t i = 0; i < num_paths; ++i) {

        if (path_clusters.path_to_cluster_index.at(i) == -1) {

            std::queue<int32_t> search_queue;
            search_queue.push(i);

            path_clusters.cluster_to_path_index.emplace_back(vector<uint32_t>());

            while (!search_queue.empty()) {

                auto cur_path = search_queue.front();

                bool is_first_visit = (path_clusters.path_to_cluster_index.at(cur_path) == -1);
                assert(is_first_visit || path_clusters.path_to_cluster_index.at(cur_path) == path_clusters.cluster_to_path_index.size() - 1);

                path_clusters.path_to_cluster_index.at(cur_path) = path_clusters.cluster_to_path_index.size() - 1;
                
                if (is_first_visit) {

                    path_clusters.cluster_to_path_index.back().emplace_back(cur_path);
                    auto connected_paths_it = connected_paths.find(cur_path);

                    if (connected_paths_it != connected_paths.end()) {

                        for (auto & path: connected_paths_it->second) {

                            if (path_clusters.path_to_cluster_index.at(path) == -1) {

                                search_queue.push(path);
                            }
                        }
                    }
                }

                search_queue.pop();
            }

            sort(path_clusters.cluster_to_path_index.back().begin(), path_clusters.cluster_to_path_index.back().end());
        }
    }

    return path_clusters;
}

double calc_read_prob(const Alignment & align, const vector<double> & quality_match_probs, const vector<double> & quality_mismatch_probs, const double indel_prob) {

    double align_path_prob = 0;

    auto & base_qualities = align.quality();
    int32_t cur_pos = 0;

    for (auto & mapping: align.path().mapping()) {

        for (auto & edit: mapping.edit()) {

            if (edit_is_match(edit)) {

                for (int32_t i = cur_pos; i < cur_pos + edit.from_length(); ++i) {

                    align_path_prob += quality_match_probs.at(int32_t(base_qualities.at(i)));
                }

            } else if (edit_is_sub(edit)) {

                for (int32_t i = cur_pos; i < cur_pos + edit.from_length(); ++i) {

                    align_path_prob += quality_mismatch_probs.at(int32_t(base_qualities.at(i)));
                }
            
            } else if (edit_is_insertion(edit)) {

                align_path_prob += edit.to_length() * indel_prob;

            } else if (edit_is_deletion(edit)) {

                align_path_prob += edit.from_length() * indel_prob;
            }
        }
    } 

    return align_path_prob;
}

AlignmentPath get_align_path(const Alignment & alignment, const gbwt::GBWT & paths_index) {
            
    auto mapping_it = alignment.path().mapping().cbegin();

    AlignmentPath align_path;
    align_path.path = paths_index.find(mapping_to_gbwt(*mapping_it));
    align_path.length = mapping_to_length(*mapping_it);
    align_path.scores.first = alignment.score();
    align_path.mapqs.first = alignment.mapping_quality();

    ++mapping_it;

    while (mapping_it != alignment.path().mapping().cend()) {

        align_path.path = paths_index.extend(align_path.path, mapping_to_gbwt(*mapping_it));
        align_path.length += mapping_to_length(*mapping_it);
        ++mapping_it;
    }

    return align_path;
}

AlignmentPath get_align_path_with_ids(const Alignment & alignment, const gbwt::GBWT & paths_index) {

    auto align_path = get_align_path(alignment, paths_index);
    align_path.path_ids = paths_index.locate(align_path.path);

    return align_path;
}

vector<gbwt::node_type> path_to_gbwt_nodes(const Path & path) {

    vector<gbwt::node_type> gbwt_nodes;
    gbwt_nodes.reserve(path.mapping_size());

    for (auto & mapping: path.mapping()) {

        gbwt_nodes.emplace_back(mapping_to_gbwt(mapping));
    }

    return gbwt_nodes;
}

void find_paired_align_paths(vector<AlignmentPath> * paired_align_paths, const AlignmentPath & start_align_path, const Alignment & end_alignment, const gbwt::GBWT & paths_index, const xg::XG & xg_index, const int32_t max_pair_distance) {

    assert(!start_align_path.path.empty());

    std::queue<AlignmentPath> paired_align_path_queue;
    paired_align_path_queue.push(start_align_path);

    assert(end_alignment.path().mapping_size() > 0);

    auto end_alignment_nodes = path_to_gbwt_nodes(end_alignment.path());

    // Perform depth-first path extension.
    while (!paired_align_path_queue.empty()) {

        auto & cur_paired_align_path = paired_align_path_queue.front();

        if (cur_paired_align_path.length > max_pair_distance) {

            paired_align_path_queue.pop();
            continue;                
        }

        auto end_alignment_nodes_it = find(end_alignment_nodes.begin(), end_alignment_nodes.end(), cur_paired_align_path.path.node);

        // Stop current extension if end node is reached.
        if (end_alignment_nodes_it != end_alignment_nodes.end()) {

            ++end_alignment_nodes_it;

            while (end_alignment_nodes_it != end_alignment_nodes.end()) {

                cur_paired_align_path.path = paths_index.extend(cur_paired_align_path.path, *end_alignment_nodes_it);
                cur_paired_align_path.length += xg_index.node_length(gbwt::Node::id(*end_alignment_nodes_it));
                ++end_alignment_nodes_it;
            }

            if (!cur_paired_align_path.path.empty() && cur_paired_align_path.length <= max_pair_distance) {

                paired_align_paths->emplace_back(cur_paired_align_path);

                paired_align_paths->back().path_ids = paths_index.locate(paired_align_paths->back().path);                                
                paired_align_paths->back().scores.second = end_alignment.score();
                paired_align_paths->back().mapqs.second = end_alignment.mapping_quality();           
            }

            paired_align_path_queue.pop();
            continue;            
        }
        
        auto out_edges = paths_index.edges(cur_paired_align_path.path.node);

        // End current extension if no outgoing edges exist.
        if (out_edges.empty()) {

            paired_align_path_queue.pop();
            continue;
        }

        auto out_edges_it = out_edges.begin(); 

        while (out_edges_it != out_edges.end()) {

            if (out_edges_it->first != gbwt::ENDMARKER) {

                auto extended_path = paths_index.extend(cur_paired_align_path.path, out_edges_it->first);

                // Add new extension to queue if not empty (path found).
                if (!extended_path.empty()) { 

                    AlignmentPath new_paired_align_path;
                    new_paired_align_path.path = extended_path;
                    new_paired_align_path.length = cur_paired_align_path.length + xg_index.node_length(gbwt::Node::id(out_edges_it->first));
                    new_paired_align_path.scores = cur_paired_align_path.scores;
                    new_paired_align_path.mapqs = cur_paired_align_path.mapqs;

                    paired_align_path_queue.push(new_paired_align_path);
                }
            }

            ++out_edges_it;
        }

        paired_align_path_queue.pop();
    }
}

vector<AlignmentPath> get_paired_align_paths(const Alignment & alignment_1, const Alignment & alignment_2, const gbwt::GBWT & paths_index, const xg::XG & xg_index, const int32_t max_pair_distance) {

    vector<AlignmentPath> paired_align_paths;

    function<size_t(const int64_t)> node_length_func = [&xg_index](const int64_t node_id) { return xg_index.node_length(node_id); };

    auto align_path_1 = get_align_path(alignment_1, paths_index);
    if (!align_path_1.path.empty()) {

        const Alignment alignment_2_rc = reverse_complement_alignment(alignment_2, node_length_func);
        find_paired_align_paths(&paired_align_paths, align_path_1, alignment_2_rc, paths_index, xg_index, max_pair_distance);
    }

    auto align_path_2 = get_align_path(alignment_2, paths_index);
    if (!align_path_2.path.empty()) {

        const Alignment alignment_1_rc = reverse_complement_alignment(alignment_1, node_length_func);
        find_paired_align_paths(&paired_align_paths, align_path_2, alignment_1_rc, paths_index, xg_index, max_pair_distance);
    }

    return paired_align_paths;
}

pair<double, vector<double> > calculate_paired_path_probs(const vector<AlignmentPath> & align_paths, const unordered_map<uint32_t, uint32_t> & clustered_path_index, const uint32_t frag_length_mean, const uint32_t frag_length_sd) {

    assert(!align_paths.empty());

    double map_prob = 1;

    if (align_paths.front().mapqs.first > 0 && align_paths.front().mapqs.second > 0) {

        map_prob = 1 - (1 - phred_to_prob(align_paths.front().mapqs.first)) * (1 - phred_to_prob(align_paths.front().mapqs.second));
    
    } else {

        return make_pair(map_prob, vector<double>(clustered_path_index.size(), 0));        
    }

    vector<double> align_paths_probs;
    align_paths_probs.reserve(align_paths.size());

    if (align_paths.size() == 1) {

        align_paths_probs.emplace_back(1 - map_prob);
    
    } else {

        double sum_scores = 0;

        for (auto & align_path: align_paths) {

            sum_scores += align_path.scores.first + align_path.scores.second;
        }

        double sum_probs = 0;

        for (auto & align_path: align_paths) {

            align_paths_probs.emplace_back((align_path.scores.first + align_path.scores.second) / sum_scores * normal_pdf<double>(align_path.length, frag_length_mean, frag_length_sd));
            sum_probs += align_paths_probs.back();
        }

        for (auto & probs: align_paths_probs) {

            probs /= sum_probs;
            probs *= (1 - map_prob);
        }
    }

    vector<double> paired_path_probs(clustered_path_index.size(), 0);

    for (size_t i = 0; i < align_paths.size(); ++i) {

        for (auto & path: align_paths.at(i).path_ids) {

            paired_path_probs.at(clustered_path_index.at(path)) = align_paths_probs.at(i) / align_paths.at(i).path_ids.size();
        }
    }

    return make_pair(map_prob, paired_path_probs);
} 

inline bool doubleCompare(const double a, const double b) {

    assert(std::isfinite(a));
    assert(std::isfinite(b));

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

bool is_path_probs_identical(const pair<double, vector<double> > & path_probs_1, const pair<double, vector<double> > & path_probs_2) {

    if (doubleCompare(path_probs_1.first, path_probs_2.first)) {

        if (path_probs_1.second.size() == path_probs_2.second.size()) {

            for (size_t i = 0; i < path_probs_1.second.size(); ++i) {

                if (!doubleCompare(path_probs_1.second.at(i), path_probs_2.second.at(i))) {

                    return false;
                }
            }

            return true;
        }
    } 

    return false;
}

void help_probs(char** argv) {
    cerr << "\nusage: " << argv[0] << " probs [options] <graph.xg> <paths.gbwt> <alignments.gam(p)> > align_path_probs.txt" << endl
         << "options:" << endl
         << "    -m, --multipath            input is multipath alignment format (GAMP)" << endl
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
    
    bool is_multipath = false;
    int32_t num_threads = 1;
    bool show_progress = false;

    int32_t c;
    optind = 2;

    while (true) {
        static struct option long_options[] =
            {
                {"multipath", no_argument, 0, 'm'},
                {"threads", no_argument, 0, 't'},
                {"progress", no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int32_t option_index = 0;
        c = getopt_long(argc, argv, "mt:ph?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'm':
            is_multipath = true;
            break;

        case 't':
            num_threads = parse<int>(optarg);
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
    omp_set_num_threads(num_threads);

    double time1 = gcsa::readTimer();

    unique_ptr<xg::XG> xg_index = vg::io::VPKG::load_one<xg::XG>(get_input_file_name(optind, argc, argv));

    double time2 = gcsa::readTimer();
    cerr << "Load XG " << time2 - time1 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    unique_ptr<gbwt::GBWT> paths_index = vg::io::VPKG::load_one<gbwt::GBWT>(get_input_file_name(optind, argc, argv));

    double time3 = gcsa::readTimer();
    cerr << "Load GBWT " << time3 - time2 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    const auto num_paths = paths_index->metadata.haplotype_count;

    vector<unordered_map<int32_t, unordered_set<int32_t> > > connected_paths_threads(num_threads);
    vector<vector<vector<AlignmentPath> > > paired_align_paths_threads(num_threads);

    get_input_file(optind, argc, argv, [&](istream& in) {

        vg::io::for_each_interleaved_pair_parallel<Alignment>(in, [&](Alignment& alignment_1, Alignment& alignment_2) {

            if (alignment_1.has_path() > 0 && alignment_2.has_path() > 0 && alignment_1.mapping_quality() > 0 && alignment_2.mapping_quality() > 0) {

                auto paired_align_paths = get_paired_align_paths(alignment_1, alignment_2, *paths_index, *xg_index, frag_length_mean + 10 * frag_length_sd);

#ifdef debug
                cerr << endl;
                cerr << pb2json(alignment_1) << endl;
                cerr << pb2json(alignment_2) << endl;

                cerr << get_align_path_with_ids(alignment_1, *paths_index) << endl;
                cerr << get_align_path_with_ids(reverse_complement_alignment(alignment_1, [&xg_index](const int64_t node_id) { return xg_index->node_length(node_id); }), *paths_index) << endl;
                cerr << get_align_path_with_ids(alignment_2, *paths_index) << endl;
                cerr << get_align_path_with_ids(reverse_complement_alignment(alignment_2, [&xg_index](const int64_t node_id) { return xg_index->node_length(node_id); }), *paths_index) << endl;
          
                for (auto & align_path: paired_align_paths) {

                    cerr << align_path << endl;
                }
#endif 

                if (!paired_align_paths.empty()) {

                    auto anchor_path_id = paired_align_paths.front().path_ids.front();

                    for (auto & align_path: paired_align_paths) {

                        for (auto & path_id: align_path.path_ids) {

                            if (anchor_path_id != path_id) {

                                connected_paths_threads.at(omp_get_thread_num())[anchor_path_id].emplace(path_id);
                                connected_paths_threads.at(omp_get_thread_num())[path_id].emplace(anchor_path_id);
                            }
                        }
                    }

                    paired_align_paths_threads.at(omp_get_thread_num()).emplace_back(paired_align_paths);

                    if (paired_align_paths_threads.at(omp_get_thread_num()).size() % 1000000 == 0) {

                        cerr << omp_get_thread_num() << ": " << paired_align_paths_threads.at(omp_get_thread_num()).size() << endl;        
                    }
                }
            }
        });
    });

    double time5 = gcsa::readTimer();
    cerr << "Found paired alignment paths " << time5 - time3 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    for (size_t i = 1; i < connected_paths_threads.size(); ++i) {

        for (auto & connected_path_clusters: connected_paths_threads.at(i)) {

            auto connected_paths_threads_it = connected_paths_threads.front().emplace(connected_path_clusters.first, unordered_set<int32_t>());
            for (auto & paths: connected_path_clusters.second) {

                connected_paths_threads_it.first->second.emplace(paths);
            }
        }
    }

    double time6 = gcsa::readTimer();
    cerr << "Merged connected path threads " << time6 - time5 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    auto path_clusters = find_path_clusters(connected_paths_threads.front(), num_paths);

    double time62 = gcsa::readTimer();
    cerr << "Found path clusters " << time62 - time6 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
 
    vector<vector<vector<AlignmentPath> > > clustered_paired_align_paths(path_clusters.cluster_to_path_index.size());

    for (auto & paired_align_paths: paired_align_paths_threads) {

        for (auto & align_path: paired_align_paths) {

            clustered_paired_align_paths.at(path_clusters.path_to_cluster_index.at(align_path.front().path_ids.front())).push_back(move(align_path));
        }
    }

    double time7 = gcsa::readTimer();
    cerr << "Clustered paired alignment paths " << time7 - time6 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    #pragma omp parallel
    { 
    
        #pragma omp for
        for (size_t i = 0; i < clustered_paired_align_paths.size(); ++i) {

            sort(clustered_paired_align_paths.at(i).begin(), clustered_paired_align_paths.at(i).end(), AlignmentPathsSorter);
        }
    }

    double time8 = gcsa::readTimer();
    cerr << "Sorted paired alignment paths " << time8 - time7 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    for (size_t i = 0; i < clustered_paired_align_paths.size(); ++i) {

        vector<vector<AlignmentPath> > & paired_align_paths = clustered_paired_align_paths.at(i);

        if (paired_align_paths.empty()) {

            continue;
        }

        unordered_map<uint32_t, uint32_t> clustered_path_index;

        cout << "#\nx 0";
        for (auto & path_id: path_clusters.cluster_to_path_index.at(i)) {

            assert(clustered_path_index.emplace(path_id, clustered_path_index.size()).second);
            cout << " " << path_id + 1;
        }
        cout << endl;

        int32_t num_paired_paths = 1;
        auto paired_path_probs = calculate_paired_path_probs(paired_align_paths.front(), clustered_path_index, frag_length_mean, frag_length_sd);

        for (size_t j = 1; j < paired_align_paths.size(); ++j) {

            auto paired_path_probs_new = calculate_paired_path_probs(paired_align_paths.at(j), clustered_path_index, frag_length_mean, frag_length_sd);

            if (is_path_probs_identical(paired_path_probs, paired_path_probs_new)) {

                num_paired_paths++;
            
            } else {

                cout << num_paired_paths << " " << paired_path_probs.first;
                for (auto & prob: paired_path_probs.second) {

                    cout << " " << prob;
                }
                cout << endl;

                num_paired_paths = 1;
                paired_path_probs = paired_path_probs_new;
            }
        }

        cout << num_paired_paths << " " << paired_path_probs.first;
        for (auto & prob: paired_path_probs.second) {

            cout << " " << prob;
        }
        cout << endl;
    }
 
    double time9 = gcsa::readTimer();
    cerr << "Wrote output " << time9 - time8 << " seconds, " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;

    return 0;
}

// Register subcommand
static Subcommand vg_probs("probs", "calculate read path probabilities", main_probs);

