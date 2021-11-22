#ifndef HALFEDGE
#define HALFEDGE

#include "../externel/obj_loader/obj_loader.h"

#include "vector3.h"
using namespace std;
#include <vector>
#include <string>
#include <list>
#include <stack>

struct HE_Edge;
struct HE_Face;

struct HE_Edge {
	HE_Edge* e_pair;  //对偶边
	HE_Edge* e_succ;  //后继边
	Vector3 startPoint;
	Vector3 endPoint;
	HE_Face* e_face;  //右侧面
	bool isUsed;
};

struct HE_Face {
	vector<Vector3> f_verts; //组成面的点
	vector<HE_Edge*> f_edges;         //组成面的边
	Vector3 normal;
	bool isGenerated;
	bool isUsed;
	int isSunkFace;
};

class HE_Mesh {
public:
	HE_Mesh();
	~HE_Mesh();
	bool LoadFromObj(objl::Loader Loader, const string& filename);
	bool IsFaceUesd();
	void ResetEdgeState();
	void ResetFaceState();
	void addSunkFace(HE_Face* f);
public:
	vector<Vector3> vertex;
	vector<HE_Edge*> edge;
	vector<HE_Face*> face;
	vector<HE_Face*> sunkFace;
};

bool IsDualEdge(HE_Edge* edge1, HE_Edge* edge2);
bool IsSunk(HE_Face* tri1, HE_Face* tri2, HE_Edge* coEdges);
void addDataToBuffer(HE_Mesh& mesh, float* preVe, float*& Ve, int size, float& newSize);

bool IsInTriangle(HE_Face* tri, Vector3& P);
HE_Mesh* createHeMesh(vector<HE_Face*>& face);
void splitCavity(HE_Mesh* mesh, vector<HE_Mesh*>& meshList);
void firstFilteringCavity(vector<HE_Mesh*>& meshList, vector<HE_Mesh*>& newMeshList);

Vector3 getMidPoint(HE_Face* face);

void drawFromeMesh(vector<HE_Mesh*>meshList, float*& cavityVe, int& size);
void getCavityPosition(vector<HE_Mesh*>& meshList, vector<Vector3>& positions);

vector<string> split(const string& str, const string& pattern);
void calculateDeepCenter();

#endif