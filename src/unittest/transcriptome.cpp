/// \file transcriptome.cpp
///  
/// unit tests for transcriptome class
///

#include <stdio.h>
#include <iostream>

#include "gbwt/dynamic_gbwt.h"
#include "bdsg/hash_graph.hpp"

#include "../transcriptome.hpp"

#include "catch.hpp"

namespace vg {
    namespace unittest {

        vector<uint64_t> handleVectorToIntVector(const vector<handle_t> & handle_vector) {

            vector<uint64_t> int_vector;
            int_vector.reserve(handle_vector.size());

            for (auto & handle: handle_vector) {

                int_vector.emplace_back(bdsg::as_integer(handle));
            }

            return int_vector;
        }

        TEST_CASE("Transcriptome can add splice-junctions and project transcripts", "[transcriptome]") {
         
            unique_ptr<MutablePathDeletableHandleGraph> graph(new bdsg::HashGraph);
               
            handle_t node1 = graph->create_handle("AAAA");
            handle_t node2 = graph->create_handle("CC");
            handle_t node3 = graph->create_handle("G");
            handle_t node4 = graph->create_handle("TTTTTTTT");
            handle_t node5 = graph->create_handle("CCC");
            handle_t node6 = graph->create_handle("AAAA");
            
            graph->create_edge(node1, node2);
            graph->create_edge(node1, node3);
            graph->create_edge(node2, node4);
            graph->create_edge(node3, node4);
            graph->create_edge(node4, node5);
            graph->create_edge(node4, node6);
            graph->create_edge(node5, node6);

            path_handle_t path1 = graph->create_path_handle("path1");
            graph->append_step(path1, node1);
            graph->append_step(path1, node2);
            graph->append_step(path1, node4);
            graph->append_step(path1, node6);

            path_handle_t path2 = graph->create_path_handle("path2");
            graph->append_step(path2, node1);
            graph->append_step(path2, node2);
            graph->append_step(path2, node4);
            graph->append_step(path2, node5);
            graph->append_step(path2, node6);

            path_handle_t path3 = graph->create_path_handle("path3");
            graph->append_step(path3, node1);
            graph->append_step(path3, node3);
            graph->append_step(path3, node4);
            graph->append_step(path3, node5);
            graph->append_step(path3, node6);

            Transcriptome transcriptome(move(graph));
            REQUIRE(graph == nullptr);

            REQUIRE(transcriptome.size() == 0);
            REQUIRE(!transcriptome.splice_graph_node_updated());

            REQUIRE(transcriptome.splice_graph().get_node_count() == 6);
            REQUIRE(transcriptome.splice_graph().get_edge_count() == 7);
            REQUIRE(transcriptome.splice_graph().get_path_count() == 3);

            REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path1")) == 4);
            REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path2")) == 5);
            REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path3")) == 5);

            unique_ptr<gbwt::GBWT> empty_haplotype_index(new gbwt::GBWT());

            std::stringstream transcript_stream;
            transcript_stream << "path1\t.\texon\t2\t7\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
            transcript_stream << "path1\t.\texon\t9\t10\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
            transcript_stream << "path1\t.\texon\t16\t18\t.\t+\t.\ttranscript_id \"transcript1\";" << endl;
            transcript_stream << "path1\t.\texon\t2\t7\t.\t-\t.\ttranscript_id \"transcript2\";" << endl;
            transcript_stream << "path1\t.\texon\t16\t18\t.\t-\t.\ttranscript_id \"transcript2\";" << endl;

            SECTION("Transcriptome can add splice-junctions") {

                transcriptome.add_transcript_splice_junctions(transcript_stream, empty_haplotype_index);

                REQUIRE(transcriptome.size() == 0);
                REQUIRE(transcriptome.splice_graph_node_updated());

                REQUIRE(transcriptome.splice_graph().get_node_count() == 11);
                REQUIRE(transcriptome.splice_graph().get_edge_count() == 15);
                REQUIRE(transcriptome.splice_graph().get_path_count() == 3);

                REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path1")) == 9);
                REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path2")) == 10);
                REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path3")) == 10);

                transcript_stream.clear();
                transcript_stream.seekg(0,ios::beg);

                SECTION("Transcriptome can project transcripts onto reference paths") {

                    transcriptome.use_reference_paths = true;

                    transcriptome.add_transcripts(transcript_stream, *empty_haplotype_index);
                    REQUIRE(transcriptome.size() == 2);

                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().front().path) == vector<uint64_t>({20, 4, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().back().path) == vector<uint64_t>({23, 9, 5, 21}));
                }

                SECTION("Transcriptome can project transcripts onto all embedded paths") {

                    transcriptome.use_all_paths = true;

                    transcriptome.add_transcripts(transcript_stream, *empty_haplotype_index);
                    REQUIRE(transcriptome.size() == 4);

                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().front().path) == vector<uint64_t>({20, 4, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(1).path) == vector<uint64_t>({20, 6, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(2).path) == vector<uint64_t>({23, 9, 5, 21}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().back().path) == vector<uint64_t>({23, 9, 7, 21}));
                }

                SECTION("Transcriptome can project transcripts onto all embedded paths and not collapse redundant paths") {

                    transcriptome.use_all_paths = true;
                    transcriptome.collapse_transcript_paths = false;

                    transcriptome.add_transcripts(transcript_stream, *empty_haplotype_index);
                    REQUIRE(transcriptome.size() == 6);

                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().front().path) == vector<uint64_t>({20, 4, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(1).path) == vector<uint64_t>({20, 4, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(2).path) == vector<uint64_t>({20, 6, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(3).path) == vector<uint64_t>({23, 9, 5, 21}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(4).path) == vector<uint64_t>({23, 9, 5, 21}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().back().path) == vector<uint64_t>({23, 9, 7, 21}));
                }

                SECTION("Transcriptome can add transcript paths to graph") {

                    transcriptome.use_all_paths = true;

                    transcriptome.add_transcripts(transcript_stream, *empty_haplotype_index);
                    REQUIRE(transcriptome.size() == 4);

                    transcriptome.embed_transcript_paths(true, true);
                    REQUIRE(transcriptome.splice_graph().get_path_count() == 7);            
                }

                SECTION("Transcriptome can remove non-transcribed nodes") {

                    transcriptome.use_all_paths = true;

                    transcriptome.add_transcripts(transcript_stream, *empty_haplotype_index);
                    REQUIRE(transcriptome.size() == 4);

                    transcriptome.remove_non_transcribed(false);

                    REQUIRE(transcriptome.splice_graph().get_node_count() == 6);
                    REQUIRE(transcriptome.splice_graph().get_edge_count() == 7);
                    REQUIRE(transcriptome.splice_graph().get_path_count() == 0);                
                }
            }

            SECTION("Transcriptome can add splice-junctions and update GBWT threads") {

                gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
                gbwt::GBWTBuilder gbwt_builder(gbwt::bit_length(gbwt::Node::encode(6, true)));

                gbwt::vector_type gbwt_thread_1(4);
                gbwt::vector_type gbwt_thread_2(5);
   
                gbwt_thread_1[0] = gbwt::Node::encode(1, false);
                gbwt_thread_1[1] = gbwt::Node::encode(2, false);
                gbwt_thread_1[2] = gbwt::Node::encode(4, false);
                gbwt_thread_1[3] = gbwt::Node::encode(6, false);

                gbwt_thread_2[0] = gbwt::Node::encode(1, false);
                gbwt_thread_2[1] = gbwt::Node::encode(2, false);
                gbwt_thread_2[2] = gbwt::Node::encode(4, false);
                gbwt_thread_2[3] = gbwt::Node::encode(5, false);
                gbwt_thread_2[4] = gbwt::Node::encode(6, false);

                gbwt::vector_type gbwt_thread_3 = gbwt_thread_1;
                gbwt_thread_3[1] = gbwt::Node::encode(3, false);

                gbwt::vector_type gbwt_thread_4 = gbwt_thread_2;
                gbwt_thread_4[1] = gbwt::Node::encode(3, false);

                gbwt_builder.insert(gbwt_thread_1, true);
                gbwt_builder.insert(gbwt_thread_2, true);
                gbwt_builder.insert(gbwt_thread_3, true);
                gbwt_builder.insert(gbwt_thread_4, true);

                gbwt_builder.finish();

                std::stringstream gbwt_stream;
                gbwt_builder.index.serialize(gbwt_stream);

                unique_ptr<gbwt::GBWT> haplotype_index(new gbwt::GBWT());
                
                haplotype_index->load(gbwt_stream);
                REQUIRE(haplotype_index->bidirectional());

                transcriptome.add_transcript_splice_junctions(transcript_stream, haplotype_index);

                REQUIRE(transcriptome.size() == 0);
                REQUIRE(transcriptome.splice_graph_node_updated());

                REQUIRE(transcriptome.splice_graph().get_node_count() == 11);
                REQUIRE(transcriptome.splice_graph().get_edge_count() == 15);
                REQUIRE(transcriptome.splice_graph().get_path_count() == 3);

                REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path1")) == 9);
                REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path2")) == 10);
                REQUIRE(transcriptome.splice_graph().get_step_count(transcriptome.splice_graph().get_path_handle("path3")) == 10);

                transcript_stream.clear();
                transcript_stream.seekg(0,ios::beg);

                SECTION("Transcriptome can project transcripts onto GBWT threads") {

                    transcriptome.add_transcripts(transcript_stream, *haplotype_index);
                    REQUIRE(transcriptome.size() == 4);

                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().front().path) == vector<uint64_t>({20, 4, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(1).path) == vector<uint64_t>({20, 6, 8, 16, 22}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().at(2).path) == vector<uint64_t>({23, 9, 5, 21}));
                    REQUIRE(handleVectorToIntVector(transcriptome.transcript_paths().back().path) == vector<uint64_t>({23, 9, 7, 21}));
                }
            }
        }
    }
}

