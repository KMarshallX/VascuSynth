/*=========================================================================
	Program: VascuSynth
	Module: Export.h

	Ground-truth export added to the VascuSynth fork: centreline graph as
	GraphML and voxel volumes (skeleton, vessel mask, vessel image, noisy
	images) as NIfTI.

	Coordinate conventions (shared by every file written here):
	  - Voxel index (i, j, k) is NIfTI array axis order; voxel centres sit at
	    integer indices.
	  - A tree position p (oxygenation-map voxels) maps to the image voxel
	    index p * mapVoxelWidth / imageVoxelWidth.
	  - World coordinates are in mm: world = voxel index * imageVoxelWidth,
	    with origin 0 and no rotation. VOXEL_WIDTH and the imageVoxelWidth
	    command-line argument are interpreted as mm.
=========================================================================*/

#ifndef _export_h
#define _export_h

#include <array>
#include <string>
#include <vector>

#include "TreeDrawer.h"
#include "VascularTree.h"

typedef std::array<int, 3> Voxel;

struct TreeEdge {
	int parent;
	int child;
};

/**
 * Everything the writers need, computed once per tree.
 *
 * The skeleton starts as every segment drawn as a 26-connected line between
 * its rounded end points. It is then cleaned into a thin, loop-free voxel
 * tree (foreground 26 / background 6 digital topology): simple points are
 * thinned away, farthest from the true centreline first; a loop is broken
 * by cutting a voxel whose neighbours split locally but stay joined
 * globally, preferring spots where unrelated branches touch; a tiny ring
 * that cannot be cut without detaching a branch is filled instead. Root and
 * terminal voxels are always kept; bifurcation voxels are kept unless a
 * loop can only be cut through one.
 */
struct GroundTruth {
	std::vector<int> nodes;                     // tree walk order from the root
	std::vector<TreeEdge> edges;                // parent -> child
	std::vector<unsigned char> skeleton;        // cleaned, 0/1, TreeDrawer layout
	std::vector<std::vector<Voxel> > edgeVoxels; // each edge's path through the cleaned skeleton
	long rawVoxels;                             // before cleaning
	long cleanVoxels;                           // after cleaning
	long tunnelsOpened;                         // spurious loops removed while cleaning
};

GroundTruth buildGroundTruth(VascularTree *vt, TreeDrawer *td);

/**
 * Writes a 3D uint8 volume as NIfTI (.nii or .nii.gz) with isotropic
 * spacing in mm and an identity-oriented affine diag(spacing, spacing,
 * spacing). `data` uses TreeDrawer's layout: data[k + dim[2]*(j + dim[1]*i)].
 */
void writeNiftiVolume(const unsigned char *data, const int dim[3], double spacingMm, const std::string &path);

/**
 * Writes the binary vessel mask: 1 where any of the voxel's eight
 * sub-samples is inside a vessel (image value > 0), plus every skeleton
 * voxel, so the mask always contains the whole skeleton.
 */
void writeVesselMaskNifti(TreeDrawer *td, const GroundTruth &gt, const std::string &path);

/** Writes the cleaned binary skeleton. */
void writeSkeletonNifti(TreeDrawer *td, const GroundTruth &gt, const std::string &path);

/**
 * Writes the tree as an undirected GraphML graph compatible with SkelHub's
 * GraphML evaluator (node `voxel_pos`, edge `centerline_voxel_points`,
 * optional world `X/Y/Z` and `centerline_world_points` in mm). Edges are
 * written parent -> child and also carry explicit `parent` / `child`
 * attributes, radius (mm and voxels), flow, length and the edge's voxels in
 * the cleaned skeleton.
 */
void writeCenterlineGraphML(VascularTree *vt, TreeDrawer *td, const GroundTruth &gt,
                            const std::string &path, const std::string &paramFile);

#endif
