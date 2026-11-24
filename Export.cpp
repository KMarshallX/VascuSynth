/*=========================================================================
	Program: VascuSynth
	Module: Export.cpp

	Ground-truth export added to the VascuSynth fork. See Export.h for the
	coordinate conventions shared by every file written here.
=========================================================================*/

#include "Export.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <queue>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>

#include "itkImage.h"
#include "itkImageFileWriter.h"

using namespace std;

namespace {

typedef array<double, 3> Point;

/** Number formatting used for every coordinate, so equal values print identically. */
string num(double value) {
	char buffer[64];
	snprintf(buffer, sizeof(buffer), "%.12g", value);
	return string(buffer);
}

string jsonPoint(const Point &p) {
	return "[" + num(p[0]) + "," + num(p[1]) + "," + num(p[2]) + "]";
}

string jsonVoxel(const Voxel &v) {
	ostringstream os;
	os << "[" << v[0] << "," << v[1] << "," << v[2] << "]";
	return os.str();
}

string xmlEscape(const string &text) {
	string out;
	for (size_t i = 0; i < text.size(); i++) {
		switch (text[i]) {
			case '&': out += "&amp;"; break;
			case '<': out += "&lt;"; break;
			case '>': out += "&gt;"; break;
			case '"': out += "&quot;"; break;
			case '\'': out += "&apos;"; break;
			default: out += text[i];
		}
	}
	return out;
}

/** Tree position (oxygenation-map voxels) -> image voxel index. */
Point toImageVoxel(VascularTree *vt, TreeDrawer *td, int node) {
	double *pos = vt->nt.getPos(node);
	double scale = td->mapVoxelWidth / td->imageVoxelWidth;
	Point p = {{pos[0] * scale, pos[1] * scale, pos[2] * scale}};
	return p;
}

Point toWorldMm(TreeDrawer *td, const Point &voxel) {
	Point p = {{voxel[0] * td->imageVoxelWidth, voxel[1] * td->imageVoxelWidth, voxel[2] * td->imageVoxelWidth}};
	return p;
}

/**
 * 26-connected line between the rounded end points, matching SkelHub's
 * `_rounded_segment`: round both ends (half to even, like numpy.rint),
 * take max|delta| steps and round each intermediate point.
 */
vector<Voxel> rasterSegment(const Point &a, const Point &b) {
	Voxel start, end, delta;
	int steps = 0;
	for (int axis = 0; axis < 3; axis++) {
		start[axis] = (int) nearbyint(a[axis]);
		end[axis] = (int) nearbyint(b[axis]);
		delta[axis] = end[axis] - start[axis];
		steps = max(steps, abs(delta[axis]));
	}
	vector<Voxel> points;
	if (steps == 0) {
		points.push_back(start);
		return points;
	}
	for (int s = 0; s <= steps; s++) {
		double t = (double) s / steps;
		Voxel v;
		for (int axis = 0; axis < 3; axis++) {
			v[axis] = (int) nearbyint(start[axis] + delta[axis] * t);
		}
		if (points.empty() || v != points.back()) {
			points.push_back(v);
		}
	}
	return points;
}

bool inImage(TreeDrawer *td, const Voxel &v) {
	for (int axis = 0; axis < 3; axis++) {
		if (v[axis] < 0 || v[axis] >= td->dim[axis]) {
			return false;
		}
	}
	return true;
}

/** Segment voxels that fall inside the image (the ones drawn into the skeleton). */
vector<Voxel> segmentVoxels(VascularTree *vt, TreeDrawer *td, const TreeEdge &e, long &clipped) {
	vector<Voxel> all = rasterSegment(toImageVoxel(vt, td, e.parent), toImageVoxel(vt, td, e.child));
	vector<Voxel> kept;
	for (size_t i = 0; i < all.size(); i++) {
		if (inImage(td, all[i])) {
			kept.push_back(all[i]);
		} else {
			clipped++;
		}
	}
	return kept;
}

int findRoot(VascularTree *vt) {
	int size = (int) vt->nt.nodes.size();
	for (int i = 0; i < size; i++) {
		if (vt->nt.getType(i) == NodeTable::ROOT) {
			return i;
		}
	}
	throw string("Unable to find the root node; the graph cannot be exported.");
}

/**
 * Walks the tree from the root exactly like the GXL writer: the root has one
 * (left) child, bifurcations have two, terminals none. Returns nodes in walk
 * order and edges parent -> child.
 */
void collectTree(VascularTree *vt, vector<int> &nodes, vector<TreeEdge> &edges) {
	vector<int> stack;
	stack.push_back(findRoot(vt));
	while (!stack.empty()) {
		int node = stack.back();
		stack.pop_back();
		nodes.push_back(node);
		int type = vt->nt.getType(node);
		if (type == NodeTable::TERM) {
			continue;
		}
		vector<int> children;
		children.push_back(vt->nt.getLeftChild(node));
		if (type != NodeTable::ROOT) {
			children.push_back(vt->nt.getRightChild(node));
		}
		for (size_t c = 0; c < children.size(); c++) {
			TreeEdge e = {node, children[c]};
			edges.push_back(e);
		}
		for (size_t c = children.size(); c-- > 0;) {
			stack.push_back(children[c]);
		}
	}
}

string nodeKind(VascularTree *vt, int node) {
	int type = vt->nt.getType(node);
	if (type == NodeTable::ROOT) return "root";
	if (type == NodeTable::TERM) return "terminal";
	if (type == NodeTable::BIF) return "bifurcation";
	return "unknown";
}

}  // namespace

void writeNiftiVolume(const unsigned char *data, const int dim[3], double spacingMm, const string &path) {
	typedef itk::Image<unsigned char, 3> VolumeType;
	VolumeType::Pointer image = VolumeType::New();

	VolumeType::SizeType size;
	VolumeType::IndexType start;
	for (int axis = 0; axis < 3; axis++) {
		size[axis] = dim[axis];
		start[axis] = 0;
	}
	VolumeType::RegionType region;
	region.SetSize(size);
	region.SetIndex(start);
	image->SetRegions(region);
	image->Allocate();

	VolumeType::SpacingType spacing;
	spacing.Fill(spacingMm);
	image->SetSpacing(spacing);

	VolumeType::PointType origin;
	origin.Fill(0.0);
	image->SetOrigin(origin);

	// ITK works in LPS and its NIfTI writer flips x and y to RAS. Pre-flipping
	// here makes the written affine diag(spacing, spacing, spacing), so the
	// world coordinate of voxel index v is v * spacing (as in the GraphML).
	VolumeType::DirectionType direction;
	direction.SetIdentity();
	direction[0][0] = -1.0;
	direction[1][1] = -1.0;
	image->SetDirection(direction);

	VolumeType::IndexType index;
	for (int i = 0; i < dim[0]; i++) {
		for (int j = 0; j < dim[1]; j++) {
			for (int k = 0; k < dim[2]; k++) {
				index[0] = i;
				index[1] = j;
				index[2] = k;
				image->SetPixel(index, data[k + dim[2] * (j + dim[1] * i)]);
			}
		}
	}

	typedef itk::ImageFileWriter<VolumeType> WriterType;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(path);
	writer->SetInput(image);
	writer->SetUseCompression(true);
	try {
		writer->Update();
	} catch (itk::ExceptionObject &excp) {
		throw string("Could not write NIfTI file ") + path + ": " + excp.GetDescription();
	}
}


namespace {

/* ---------- local digital topology (foreground 26, background 6) ---------- */

/** 3x3x3 neighbourhood cell index; the centre is 13. */
inline int cell(int dx, int dy, int dz) {
	return (dx + 1) * 9 + (dy + 1) * 3 + (dz + 1);
}

struct CubeTables {
	vector<int> adj26[27];  // 26-neighbours inside the cube, centre excluded
	vector<int> adj6[27];   // 6-neighbours inside the cube, centre excluded
	bool inN18[27];         // cells sharing a face or an edge with the centre
	bool inN6[27];          // cells sharing a face with the centre
	CubeTables() {
		for (int a = 0; a < 27; a++) {
			int ax = a / 9 - 1, ay = (a / 3) % 3 - 1, az = a % 3 - 1;
			int nonzero = (ax != 0) + (ay != 0) + (az != 0);
			inN18[a] = (a != 13) && nonzero <= 2;
			inN6[a] = (a != 13) && nonzero == 1;
			for (int b = 0; b < 27; b++) {
				if (a == b || a == 13 || b == 13) continue;
				int bx = b / 9 - 1, by = (b / 3) % 3 - 1, bz = b % 3 - 1;
				int dx = abs(ax - bx), dy = abs(ay - by), dz = abs(az - bz);
				if (max(dx, max(dy, dz)) == 1) adj26[a].push_back(b);
				if (dx + dy + dz == 1) adj6[a].push_back(b);
			}
		}
	}
};

const CubeTables &tables() {
	static CubeTables t;
	return t;
}

/** Number of 26-connected foreground components among the 26 neighbours. */
int topoT26(const bool fg[27]) {
	const CubeTables &t = tables();
	bool seen[27] = {false};
	int components = 0;
	for (int start = 0; start < 27; start++) {
		if (start == 13 || !fg[start] || seen[start]) continue;
		components++;
		vector<int> stack(1, start);
		seen[start] = true;
		while (!stack.empty()) {
			int c = stack.back();
			stack.pop_back();
			for (size_t n = 0; n < t.adj26[c].size(); n++) {
				int b = t.adj26[c][n];
				if (fg[b] && !seen[b]) {
					seen[b] = true;
					stack.push_back(b);
				}
			}
		}
	}
	return components;
}

/** Number of 6-connected background components in N18 that touch a face of the centre. */
int topoT6bar(const bool fg[27]) {
	const CubeTables &t = tables();
	bool seen[27] = {false};
	int components = 0;
	for (int start = 0; start < 27; start++) {
		if (!t.inN6[start] || fg[start] || seen[start]) continue;
		components++;
		vector<int> stack(1, start);
		seen[start] = true;
		while (!stack.empty()) {
			int c = stack.back();
			stack.pop_back();
			for (size_t n = 0; n < t.adj6[c].size(); n++) {
				int b = t.adj6[c][n];
				if (t.inN18[b] && !fg[b] && !seen[b]) {
					seen[b] = true;
					stack.push_back(b);
				}
			}
		}
	}
	return components;
}

struct Volume {
	int d[3];
	vector<unsigned char> v;
	Volume(const int dim[3]) {
		d[0] = dim[0]; d[1] = dim[1]; d[2] = dim[2];
		v.assign((size_t) d[0] * d[1] * d[2], 0);
	}
	long index(int i, int j, int k) const { return k + (long) d[2] * (j + (long) d[1] * i); }
	long index(const Voxel &x) const { return index(x[0], x[1], x[2]); }
	Voxel voxel(long n) const {
		Voxel x = {{(int) (n / ((long) d[1] * d[2])), (int) ((n / d[2]) % d[1]), (int) (n % d[2])}};
		return x;
	}
	bool inside(int i, int j, int k) const {
		return i >= 0 && j >= 0 && k >= 0 && i < d[0] && j < d[1] && k < d[2];
	}
	bool get(int i, int j, int k) const { return inside(i, j, k) && v[index(i, j, k)]; }
	void cube(long n, bool fg[27]) const {
		Voxel x = voxel(n);
		for (int dx = -1; dx <= 1; dx++)
			for (int dy = -1; dy <= 1; dy++)
				for (int dz = -1; dz <= 1; dz++)
					fg[cell(dx, dy, dz)] = get(x[0] + dx, x[1] + dy, x[2] + dz);
	}
};

double pointSegmentDistance(const Point &p, const Point &a, const Point &b) {
	double ab[3], ap[3], len2 = 0, t = 0;
	for (int i = 0; i < 3; i++) {
		ab[i] = b[i] - a[i];
		ap[i] = p[i] - a[i];
		len2 += ab[i] * ab[i];
		t += ab[i] * ap[i];
	}
	t = len2 > 0 ? max(0.0, min(1.0, t / len2)) : 0.0;
	double d2 = 0;
	for (int i = 0; i < 3; i++) {
		double diff = ap[i] - t * ab[i];
		d2 += diff * diff;
	}
	return sqrt(d2);
}

Voxel roundedNodeVoxel(VascularTree *vt, TreeDrawer *td, int node) {
	Point p = toImageVoxel(vt, td, node);
	Voxel x = {{(int) nearbyint(p[0]), (int) nearbyint(p[1]), (int) nearbyint(p[2])}};
	return x;
}

/** Shortest 26-connected path (Euclidean step lengths) between two skeleton voxels. */
bool skeletonPath(const Volume &s, const Voxel &from, const Voxel &to, vector<Voxel> &path) {
	long src = s.index(from), dst = s.index(to);
	if (!s.v[src] || !s.v[dst]) return false;
	typedef pair<double, long> Item;
	priority_queue<Item, vector<Item>, greater<Item> > queue;
	unordered_map<long, double> dist;
	unordered_map<long, long> prev;
	dist[src] = 0;
	queue.push(Item(0, src));
	while (!queue.empty()) {
		Item top = queue.top();
		queue.pop();
		if (top.second == dst) break;
		if (top.first > dist[top.second]) continue;
		Voxel x = s.voxel(top.second);
		for (int dx = -1; dx <= 1; dx++)
			for (int dy = -1; dy <= 1; dy++)
				for (int dz = -1; dz <= 1; dz++) {
					if (!dx && !dy && !dz) continue;
					if (!s.get(x[0] + dx, x[1] + dy, x[2] + dz)) continue;
					long n = s.index(x[0] + dx, x[1] + dy, x[2] + dz);
					double nd = top.first + sqrt((double) (dx * dx + dy * dy + dz * dz));
					unordered_map<long, double>::iterator it = dist.find(n);
					if (it == dist.end() || nd < it->second) {
						dist[n] = nd;
						prev[n] = top.second;
						queue.push(Item(nd, n));
					}
				}
	}
	if (!dist.count(dst)) return false;
	path.clear();
	for (long n = dst;; n = prev[n]) {
		path.push_back(s.voxel(n));
		if (n == src) break;
	}
	reverse(path.begin(), path.end());
	return true;
}

/** The voxel itself if it is in the skeleton, else the nearest skeleton voxel within 2 voxels. */
Voxel nearestSkeletonVoxel(const Volume &s, const Voxel &x) {
	if (s.get(x[0], x[1], x[2])) return x;
	Voxel best = x;
	int bestD = 1 << 30;
	for (int dx = -2; dx <= 2; dx++)
		for (int dy = -2; dy <= 2; dy++)
			for (int dz = -2; dz <= 2; dz++) {
				int d = dx * dx + dy * dy + dz * dz;
				if (d < bestD && s.get(x[0] + dx, x[1] + dy, x[2] + dz)) {
					bestD = d;
					Voxel y = {{x[0] + dx, x[1] + dy, x[2] + dz}};
					best = y;
				}
			}
	return best;
}

/** Skeleton voxels whose removal disconnects the 26-adjacency graph (Tarjan, iterative). */
vector<char> articulationPoints(const Volume &s, const vector<long> &voxels) {
	unordered_map<long, int> id;
	for (size_t i = 0; i < voxels.size(); i++) id[voxels[i]] = (int) i;
	int n = (int) voxels.size();
	vector<vector<int> > adj(n);
	for (int i = 0; i < n; i++) {
		Voxel x = s.voxel(voxels[i]);
		for (int dx = -1; dx <= 1; dx++)
			for (int dy = -1; dy <= 1; dy++)
				for (int dz = -1; dz <= 1; dz++) {
					if ((!dx && !dy && !dz) || !s.get(x[0] + dx, x[1] + dy, x[2] + dz)) continue;
					adj[i].push_back(id[s.index(x[0] + dx, x[1] + dy, x[2] + dz)]);
				}
	}
	vector<int> disc(n, -1), low(n, 0), parent(n, -1), next(n, 0), children(n, 0);
	vector<char> cut(n, 0);
	int timer = 0;
	for (int root = 0; root < n; root++) {
		if (disc[root] >= 0) continue;
		vector<int> stack(1, root);
		disc[root] = low[root] = timer++;
		while (!stack.empty()) {
			int u = stack.back();
			if (next[u] < (int) adj[u].size()) {
				int w = adj[u][next[u]++];
				if (disc[w] < 0) {
					parent[w] = u;
					children[u]++;
					disc[w] = low[w] = timer++;
					stack.push_back(w);
				} else if (w != parent[u]) {
					low[u] = min(low[u], disc[w]);
				}
			} else {
				stack.pop_back();
				int p = parent[u];
				if (p >= 0) {
					low[p] = min(low[p], low[u]);
					if (parent[p] >= 0 && low[u] >= disc[p]) cut[p] = 1;
				}
			}
		}
		if (children[root] > 1) cut[root] = 1;
	}
	return cut;
}

}  // namespace

GroundTruth buildGroundTruth(VascularTree *vt, TreeDrawer *td) {
	GroundTruth gt;
	collectTree(vt, gt.nodes, gt.edges);

	// 1. raw skeleton: every segment as a 26-connected line between rounded ends
	Volume s(td->dim);
	long clipped = 0;
	vector<Point> segA, segB;
	unordered_map<long, vector<int> > drawnBy;  // raw voxel -> edges that drew it
	for (size_t e = 0; e < gt.edges.size(); e++) {
		vector<Voxel> voxels = segmentVoxels(vt, td, gt.edges[e], clipped);
		for (size_t v = 0; v < voxels.size(); v++) {
			long n = s.index(voxels[v]);
			s.v[n] = 1;
			drawnBy[n].push_back((int) e);
		}
		segA.push_back(toImageVoxel(vt, td, gt.edges[e].parent));
		segB.push_back(toImageVoxel(vt, td, gt.edges[e].child));
	}
	if (clipped > 0) {
		cerr << "WARNING: " << clipped << " skeleton voxels fell outside the image and were dropped." << endl;
	}

	// 2. node voxels: root and terminals (2) are never removed; bifurcations (1)
	//    are kept while thinning but may be cut as a last resort to open a
	//    tunnel inside a tight cluster of nodes
	vector<char> anchor(s.v.size(), 0);
	for (size_t n = 0; n < gt.nodes.size(); n++) {
		Voxel x = roundedNodeVoxel(vt, td, gt.nodes[n]);
		if (!s.inside(x[0], x[1], x[2])) continue;
		char level = vt->nt.getType(gt.nodes[n]) == NodeTable::BIF ? 1 : 2;
		anchor[s.index(x)] = max(anchor[s.index(x)], level);
	}

	// 3. removal order: farthest from the true centreline first
	vector<long> order;
	vector<double> distance;
	for (long n = 0; n < (long) s.v.size(); n++) {
		if (!s.v[n]) continue;
		Voxel x = s.voxel(n);
		Point p = {{(double) x[0], (double) x[1], (double) x[2]}};
		double best = 1e300;
		for (size_t e = 0; e < segA.size(); e++) best = min(best, pointSegmentDistance(p, segA[e], segB[e]));
		order.push_back(n);
		distance.push_back(best);
	}
	gt.rawVoxels = (long) order.size();
	vector<size_t> rank(order.size());
	for (size_t i = 0; i < rank.size(); i++) rank[i] = i;
	stable_sort(rank.begin(), rank.end(), [&](size_t a, size_t b) { return distance[a] > distance[b]; });

	// 4a. thinning: remove simple points (T26 == 1 and T6bar == 1), which never
	//     changes the topology, until none is left
	gt.tunnelsOpened = 0;
	bool fg[27];
	auto thin = [&]() {
		long removed;
		do {
			removed = 0;
			for (size_t r = 0; r < rank.size(); r++) {
				long n = order[rank[r]];
				if (!s.v[n] || anchor[n]) continue;
				s.cube(n, fg);
				if (topoT26(fg) != 1 || topoT6bar(fg) != 1) continue;
				s.v[n] = 0;
				removed++;
			}
		} while (removed > 0);
	};

	// 4b. a voxel is a "contact" when it touches a voxel drawn only by other
	//     edges that do not meet it at a shared node within 2 voxels, i.e. two
	//     branches touching where the true tree has no junction
	auto nodeVoxelPoint = [&](int node) { return toImageVoxel(vt, td, node); };
	auto isContact = [&](long n) {
		const vector<int> &mine = drawnBy[n];
		Voxel x = s.voxel(n);
		Point p = {{(double) x[0], (double) x[1], (double) x[2]}};
		for (int dx = -1; dx <= 1; dx++)
			for (int dy = -1; dy <= 1; dy++)
				for (int dz = -1; dz <= 1; dz++) {
					if ((!dx && !dy && !dz) || !s.get(x[0] + dx, x[1] + dy, x[2] + dz)) continue;
					const vector<int> &theirs = drawnBy[s.index(x[0] + dx, x[1] + dy, x[2] + dz)];
					bool related = false;
					for (size_t a = 0; a < mine.size() && !related; a++)
						for (size_t b = 0; b < theirs.size() && !related; b++) {
							const TreeEdge &e = gt.edges[mine[a]], &f = gt.edges[theirs[b]];
							if (mine[a] == theirs[b]) { related = true; break; }
							int shared[4] = {e.parent, e.child, f.parent, f.child};
							for (int i = 0; i < 2 && !related; i++)
								for (int j = 2; j < 4 && !related; j++)
									if (shared[i] == shared[j]) {
										Point q = nodeVoxelPoint(shared[i]);
										double d = sqrt(pow(p[0] - q[0], 2) + pow(p[1] - q[1], 2) + pow(p[2] - q[2], 2));
										if (d <= 2.0) related = true;
									}
						}
					if (!related) return true;
				}
		return false;
	};

	// 4c. alternate: thin, then break one remaining loop by cutting a voxel
	//     whose neighbours split locally but stay joined globally, preferring
	//     contacts and then the voxel farthest from the centreline
	// every cut or fill removes at least one loop, so this always ends; the cap
	// is only a guard against an unforeseen configuration
	for (long step = 0;; step++) {
		if (step > gt.rawVoxels) {
			cerr << "WARNING: skeleton cleaning stopped early; some loops may remain." << endl;
			break;
		}
		thin();
		vector<long> voxels;
		vector<double> dist;
		for (size_t r = 0; r < rank.size(); r++) {
			long n = order[rank[r]];
			if (s.v[n]) { voxels.push_back(n); dist.push_back(distance[rank[r]]); }
		}
		vector<char> cut = articulationPoints(s, voxels);
		// priority: plain voxels before bifurcation voxels, then contacts, then distance
		long best = -1;
		int bestLevel = 3;
		bool bestContact = false;
		double bestDist = -1;
		for (size_t i = 0; i < voxels.size(); i++) {
			int level = anchor[voxels[i]];
			if (cut[i] || level == 2) continue;
			// removal must break a loop: neighbours split locally (T26 >= 2) but
			// stay joined globally (not an articulation point)
			s.cube(voxels[i], fg);
			if (topoT26(fg) < 2) continue;
			bool contact = isContact(voxels[i]);
			bool better = best < 0 || level < bestLevel ||
			              (level == bestLevel && ((contact && !bestContact) ||
			                                      (contact == bestContact && dist[i] > bestDist)));
			if (better) {
				best = voxels[i];
				bestLevel = level;
				bestContact = contact;
				bestDist = dist[i];
			}
		}
		if (best < 0) {
			// 4d. a tiny ring whose voxels all carry a branch has no voxel that
			//     can be cut without detaching something; fill it instead: add a
			//     background voxel that touches one skeleton component (T26 == 1)
			//     and separates the local background (T6bar >= 2), nearest to the
			//     centreline first. Thinning then pares the filled clump back.
			long fill = -1;
			double fillDist = 1e300;
			for (size_t i = 0; i < voxels.size(); i++) {
				Voxel x = s.voxel(voxels[i]);
				for (int dx = -1; dx <= 1; dx++)
					for (int dy = -1; dy <= 1; dy++)
						for (int dz = -1; dz <= 1; dz++) {
							int qi = x[0] + dx, qj = x[1] + dy, qk = x[2] + dz;
							if (!s.inside(qi, qj, qk) || s.get(qi, qj, qk)) continue;
							long q = s.index(qi, qj, qk);
							s.cube(q, fg);
							if (topoT26(fg) != 1 || topoT6bar(fg) < 2) continue;
							Point pq = {{(double) qi, (double) qj, (double) qk}};
							double dq = 1e300;
							for (size_t e = 0; e < segA.size(); e++) dq = min(dq, pointSegmentDistance(pq, segA[e], segB[e]));
							if (dq < fillDist) { fillDist = dq; fill = q; }
						}
			}
			if (fill < 0) break;
			s.cube(fill, fg);
			gt.tunnelsOpened += topoT6bar(fg) - 1;
			s.v[fill] = 1;
			order.push_back(fill);
			distance.push_back(fillDist);
			rank.push_back(rank.size());
			stable_sort(rank.begin(), rank.end(), [&](size_t a, size_t b) { return distance[a] > distance[b]; });
			continue;
		}
		s.cube(best, fg);
		gt.tunnelsOpened += topoT26(fg) - 1;
		s.v[best] = 0;
	}

	gt.cleanVoxels = 0;
	for (size_t n = 0; n < s.v.size(); n++) gt.cleanVoxels += s.v[n];

	// 5. each edge's voxels: its path through the cleaned skeleton
	for (size_t e = 0; e < gt.edges.size(); e++) {
		vector<Voxel> path;
		Voxel a = nearestSkeletonVoxel(s, roundedNodeVoxel(vt, td, gt.edges[e].parent));
		Voxel b = nearestSkeletonVoxel(s, roundedNodeVoxel(vt, td, gt.edges[e].child));
		bool ok = s.inside(a[0], a[1], a[2]) && s.inside(b[0], b[1], b[2]) && skeletonPath(s, a, b, path);
		if (!ok) {
			cerr << "WARNING: no skeleton path for segment n" << gt.edges[e].parent << " -> n"
			     << gt.edges[e].child << "; using its raw line." << endl;
			long ignored = 0;
			path = segmentVoxels(vt, td, gt.edges[e], ignored);
		}
		gt.edgeVoxels.push_back(path);
	}

	gt.skeleton.swap(s.v);
	return gt;
}

void writeVesselMaskNifti(TreeDrawer *td, const GroundTruth &gt, const string &path) {
	long total = (long) td->dim[0] * td->dim[1] * td->dim[2];
	vector<unsigned char> mask(total, 0);
	for (long n = 0; n < total; n++) {
		// value = 32 * (inside sub-samples) - 1, so any sub-sample inside <=> value > 0
		mask[n] = (td->image[n] > 0 || gt.skeleton[n]) ? 1 : 0;
	}
	writeNiftiVolume(&mask[0], td->dim, td->imageVoxelWidth, path);
}

void writeSkeletonNifti(TreeDrawer *td, const GroundTruth &gt, const string &path) {
	writeNiftiVolume(&gt.skeleton[0], td->dim, td->imageVoxelWidth, path);
}

void writeCenterlineGraphML(VascularTree *vt, TreeDrawer *td, const GroundTruth &gt,
                            const string &path, const string &paramFile) {
	const vector<int> &nodes = gt.nodes;
	const vector<TreeEdge> &edges = gt.edges;

	int terminals = 0;
	for (size_t n = 0; n < nodes.size(); n++) {
		if (vt->nt.getType(nodes[n]) == NodeTable::TERM) {
			terminals++;
		}
	}

	ofstream os(path.c_str());
	if (!os.is_open()) {
		throw string("Could not open GraphML file for writing: ") + path;
	}

	os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
	os << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\""
	   << " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\""
	   << " xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns"
	   << " http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n";

	// graph-level metadata
	os << "  <key id=\"g_generator\" for=\"graph\" attr.name=\"generator\" attr.type=\"string\"/>\n";
	os << "  <key id=\"g_param_file\" for=\"graph\" attr.name=\"param_file\" attr.type=\"string\"/>\n";
	os << "  <key id=\"g_spatial_unit\" for=\"graph\" attr.name=\"spatial_unit\" attr.type=\"string\"/>\n";
	os << "  <key id=\"g_map_voxel_width_mm\" for=\"graph\" attr.name=\"map_voxel_width_mm\" attr.type=\"double\"/>\n";
	os << "  <key id=\"g_image_voxel_width_mm\" for=\"graph\" attr.name=\"image_voxel_width_mm\" attr.type=\"double\"/>\n";
	os << "  <key id=\"g_image_shape\" for=\"graph\" attr.name=\"image_shape\" attr.type=\"string\"/>\n";
	os << "  <key id=\"g_coordinate_frame\" for=\"graph\" attr.name=\"coordinate_frame\" attr.type=\"string\"/>\n";
	os << "  <key id=\"g_num_terminals\" for=\"graph\" attr.name=\"num_terminals\" attr.type=\"int\"/>\n";
	// node attributes
	os << "  <key id=\"v_voxel_pos\" for=\"node\" attr.name=\"voxel_pos\" attr.type=\"string\"/>\n";
	os << "  <key id=\"v_X\" for=\"node\" attr.name=\"X\" attr.type=\"double\"/>\n";
	os << "  <key id=\"v_Y\" for=\"node\" attr.name=\"Y\" attr.type=\"double\"/>\n";
	os << "  <key id=\"v_Z\" for=\"node\" attr.name=\"Z\" attr.type=\"double\"/>\n";
	os << "  <key id=\"v_map_pos\" for=\"node\" attr.name=\"map_pos\" attr.type=\"string\"/>\n";
	os << "  <key id=\"v_kind\" for=\"node\" attr.name=\"kind\" attr.type=\"string\"/>\n";
	// edge attributes
	os << "  <key id=\"e_parent\" for=\"edge\" attr.name=\"parent\" attr.type=\"string\"/>\n";
	os << "  <key id=\"e_child\" for=\"edge\" attr.name=\"child\" attr.type=\"string\"/>\n";
	os << "  <key id=\"e_radius_mm\" for=\"edge\" attr.name=\"radius_mm\" attr.type=\"double\"/>\n";
	os << "  <key id=\"e_radius_voxels\" for=\"edge\" attr.name=\"radius_voxels\" attr.type=\"double\"/>\n";
	os << "  <key id=\"e_flow\" for=\"edge\" attr.name=\"flow\" attr.type=\"double\"/>\n";
	os << "  <key id=\"e_length_mm\" for=\"edge\" attr.name=\"length_mm\" attr.type=\"double\"/>\n";
	os << "  <key id=\"e_centerline_voxel_points\" for=\"edge\" attr.name=\"centerline_voxel_points\" attr.type=\"string\"/>\n";
	os << "  <key id=\"e_centerline_world_points\" for=\"edge\" attr.name=\"centerline_world_points\" attr.type=\"string\"/>\n";
	os << "  <key id=\"e_centerline_voxels\" for=\"edge\" attr.name=\"centerline_voxels\" attr.type=\"string\"/>\n";
	os << "  <key id=\"e_num_centerline_voxels\" for=\"edge\" attr.name=\"num_centerline_voxels\" attr.type=\"int\"/>\n";

	os << "  <graph id=\"G\" edgedefault=\"undirected\">\n";
	os << "    <data key=\"g_generator\">VascuSynth</data>\n";
	os << "    <data key=\"g_param_file\">" << xmlEscape(paramFile) << "</data>\n";
	os << "    <data key=\"g_spatial_unit\">mm</data>\n";
	os << "    <data key=\"g_map_voxel_width_mm\">" << num(td->mapVoxelWidth) << "</data>\n";
	os << "    <data key=\"g_image_voxel_width_mm\">" << num(td->imageVoxelWidth) << "</data>\n";
	os << "    <data key=\"g_image_shape\">[" << td->dim[0] << "," << td->dim[1] << "," << td->dim[2] << "]</data>\n";
	os << "    <data key=\"g_coordinate_frame\">voxel_pos and centerline_voxel_points are (i,j,k) voxel indices"
	   << " of the NIfTI volumes written alongside this file; X/Y/Z and centerline_world_points are"
	   << " index * image_voxel_width_mm</data>\n";
	os << "    <data key=\"g_num_terminals\">" << terminals << "</data>\n";

	for (size_t n = 0; n < nodes.size(); n++) {
		int node = nodes[n];
		Point voxel = toImageVoxel(vt, td, node);
		Point world = toWorldMm(td, voxel);
		double *pos = vt->nt.getPos(node);
		Point mapPos = {{pos[0], pos[1], pos[2]}};
		os << "    <node id=\"n" << node << "\">\n";
		os << "      <data key=\"v_voxel_pos\">" << jsonPoint(voxel) << "</data>\n";
		os << "      <data key=\"v_X\">" << num(world[0]) << "</data>\n";
		os << "      <data key=\"v_Y\">" << num(world[1]) << "</data>\n";
		os << "      <data key=\"v_Z\">" << num(world[2]) << "</data>\n";
		os << "      <data key=\"v_map_pos\">" << jsonPoint(mapPos) << "</data>\n";
		os << "      <data key=\"v_kind\">" << nodeKind(vt, node) << "</data>\n";
		os << "    </node>\n";
	}

	for (size_t e = 0; e < edges.size(); e++) {
		const TreeEdge &edge = edges[e];
		Point a = toImageVoxel(vt, td, edge.parent);
		Point b = toImageVoxel(vt, td, edge.child);
		const vector<Voxel> &voxels = gt.edgeVoxels[e];

		double lengthVoxels = sqrt(pow(b[0] - a[0], 2) + pow(b[1] - a[1], 2) + pow(b[2] - a[2], 2));
		if (lengthVoxels == 0.0) {
			cerr << "WARNING: segment n" << edge.parent << " -> n" << edge.child << " has zero length." << endl;
		}
		double radiusMm = vt->nt.getRadius(edge.child);

		string voxelList = "[";
		for (size_t v = 0; v < voxels.size(); v++) {
			voxelList += (v ? "," : "") + jsonVoxel(voxels[v]);
		}
		voxelList += "]";

		os << "    <edge id=\"e" << edge.child << "\" source=\"n" << edge.parent << "\" target=\"n" << edge.child << "\">\n";
		os << "      <data key=\"e_parent\">n" << edge.parent << "</data>\n";
		os << "      <data key=\"e_child\">n" << edge.child << "</data>\n";
		os << "      <data key=\"e_radius_mm\">" << num(radiusMm) << "</data>\n";
		os << "      <data key=\"e_radius_voxels\">" << num(radiusMm / td->imageVoxelWidth) << "</data>\n";
		os << "      <data key=\"e_flow\">" << num(vt->nt.getFlow(edge.child)) << "</data>\n";
		os << "      <data key=\"e_length_mm\">" << num(lengthVoxels * td->imageVoxelWidth) << "</data>\n";
		os << "      <data key=\"e_centerline_voxel_points\">[" << jsonPoint(a) << "," << jsonPoint(b) << "]</data>\n";
		os << "      <data key=\"e_centerline_world_points\">[" << jsonPoint(toWorldMm(td, a)) << ","
		   << jsonPoint(toWorldMm(td, b)) << "]</data>\n";
		os << "      <data key=\"e_centerline_voxels\">" << voxelList << "</data>\n";
		os << "      <data key=\"e_num_centerline_voxels\">" << voxels.size() << "</data>\n";
		os << "    </edge>\n";
	}

	os << "  </graph>\n";
	os << "</graphml>\n";
	os.close();
}
