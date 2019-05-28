#ifndef VG_IDENTITY_OVERLAY_HPP_INCLUDED
#define VG_IDENTITY_OVERLAY_HPP_INCLUDED

/** \file
 * identity_overlay.hpp: defines an overlay that does not modify the underlying graph
 */

#include "handle.hpp"

namespace vg {

using namespace std;

    /**
     * A HandleGraph overlay that does not modify the underlying graph but is still
     * useful for providing a uniform interface when choosing between overlays.
     */
    class IdentityOverlay : public ExpandingOverlayGraph {
    public:
        
        /// Initialize as the reverse version of another graph, optionally also
        /// complementing
        IdentityOverlay(const HandleGraph* graph);
        
        /// Default constructor -- not actually functional
        IdentityOverlay() = default;
        
        /// Default destructor
        ~IdentityOverlay() = default;
        
        //////////////////////////
        /// HandleGraph interface
        //////////////////////////
        
        // Method to check if a node exists by ID
        virtual bool has_node(id_t node_id) const;
        
        /// Look up the handle for the node with the given ID in the given orientation
        virtual handle_t get_handle(const id_t& node_id, bool is_reverse = false) const;
        
        /// Get the ID from a handle
        virtual id_t get_id(const handle_t& handle) const;
        
        /// Get the orientation of a handle
        virtual bool get_is_reverse(const handle_t& handle) const;
        
        /// Invert the orientation of a handle (potentially without getting its ID)
        virtual handle_t flip(const handle_t& handle) const;
        
        /// Get the length of a node
        virtual size_t get_length(const handle_t& handle) const;
        
        /// Get the sequence of a node, presented in the handle's local forward
        /// orientation.
        virtual string get_sequence(const handle_t& handle) const;
        
        /// Loop over all the handles to next/previous (right/left) nodes. Passes
        /// them to a callback which returns false to stop iterating and true to
        /// continue. Returns true if we finished and false if we stopped early.
        virtual bool follow_edges(const handle_t& handle, bool go_left, const function<bool(const handle_t&)>& iteratee) const;
        
        /// Loop over all the nodes in the graph in their local forward
        /// orientations, in their internal stored order. Stop if the iteratee
        /// returns false. Can be told to run in parallel, in which case stopping
        /// after a false return value is on a best-effort basis and iteration
        /// order is not defined.
        virtual void for_each_handle(const function<bool(const handle_t&)>& iteratee, bool parallel = false) const;
        
        /// Return the number of nodes in the graph
        /// TODO: can't be node_count because XG has a field named node_count.
        virtual size_t get_node_count() const;
        
        /// Return the smallest ID in the graph, or some smaller number if the
        /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t min_node_id() const;
        
        /// Return the largest ID in the graph, or some larger number if the
        /// largest ID is unavailable. Return value is unspecified if the graph is empty.
        virtual id_t max_node_id() const;
        
        ///////////////////////////////////
        /// ExpandingOverlayGraph interface
        ///////////////////////////////////
        
        /**
         * Returns the handle in the underlying graph that corresponds to a handle in the
         * overlay
         */
        virtual handle_t get_underlying_handle(const handle_t& handle) const;
        
    private:
        /// The underlying graph
        const HandleGraph* graph = nullptr;
    };
}

#endif
