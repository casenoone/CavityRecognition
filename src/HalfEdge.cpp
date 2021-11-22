#include "HalfEdge.h"

HE_Mesh::HE_Mesh() {

}

HE_Mesh::~HE_Mesh() {

}

Vector3 getMidPoint(HE_Face* face) {
	auto A = face->f_verts[0];
	auto B = face->f_verts[1];
	auto C = face->f_verts[2];

	auto x = (A.x + B.x + C.x) / 3;
	auto y = (A.y + B.y + C.y) / 3;
	auto z = (A.z + B.z + C.z) / 3;

	Vector3 newP(x, y, z);
	return newP;
}

bool IsInTriangle(HE_Face* tri, Vector3& P) {
	auto A = tri->f_verts[0];
	auto B = tri->f_verts[1];
	auto C = tri->f_verts[2];

	Vector3 v0 = C - A;
	Vector3 v1 = B - A;
	Vector3 v2 = P - A;

	float dot00 = v0.dot(v0);
	float dot01 = v0.dot(v1);
	float dot02 = v0.dot(v2);
	float dot11 = v1.dot(v1);
	float dot12 = v1.dot(v2);

	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	if (u < 0 || u > 1)
	{
		return false;
	}

	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno;
	if (v < 0 || v > 1)
	{
		return false;
	}

	return u + v <= 1;

}

bool IsSunk(HE_Face* tri1, HE_Face* tri2, HE_Edge* coEdges) {
	auto N1 = tri1->normal;
	auto N2 = tri2->normal;

	auto coEdge = coEdges->endPoint - coEdges->startPoint;
	coEdge.normalize();

	auto tangentV1 = (N1.cross(coEdge)).normalize();
	auto tangentV2 = (N2.cross(coEdge).normalize());
	auto midPoint = (coEdges->endPoint + coEdges->startPoint) / 2.0f;

	auto point1 = Vector3();
	auto point2 = Vector3();
	auto scale = 0.0001f;

	point1 = tangentV1 * scale + midPoint;
	point2 = tangentV2 * scale + midPoint;

	auto move1 = tangentV1 * scale;
	auto move2 = tangentV2 * scale;
	if (IsInTriangle(tri1, point1) == false) {
		point1 = midPoint - move1;
	}

	if (IsInTriangle(tri2, point2) == false) {
		point2 = midPoint - move2;
	}


	Ray ray1(point1, N1, 1.0f);
	Ray ray2(point2, N2, 1.0f);

	if (ray1.IsInsert(ray2) == true) {
		return true;
	}

	return false;

}

bool IsDualEdge(HE_Edge* edge1, HE_Edge* edge2) {
	if (((edge1->startPoint == edge2->startPoint) && (edge1->endPoint == edge2->endPoint)) ||
		((edge1->startPoint == edge2->endPoint) && (edge1->endPoint == edge2->startPoint))) {
		return true;
	}
	return false;
}


bool HE_Mesh::LoadFromObj(objl::Loader Loader, const string& filename) {
	bool loadout = Loader.LoadFile(filename);
	if (!loadout) {
		cout << "faild";
		return false;
	}

	list<HE_Edge*> dualList;
	for (int i = 0; i < Loader.LoadedMeshes.size(); i++) {
		objl::Mesh curMesh = Loader.LoadedMeshes[i];
		cout << "正在生成半边数据结构" << endl;


		for (int j = 0; j < curMesh.Indices.size(); j += 3)
		{
			int a, b, c;
			a = curMesh.Indices[j];
			b = curMesh.Indices[j + 1];
			c = curMesh.Indices[j + 2];

			auto n1 = Loader.LoadedVertices[a].Normal;
			auto n2 = Loader.LoadedVertices[b].Normal;
			auto n3 = Loader.LoadedVertices[c].Normal;

			Vector3 N1(n1.X, n1.Y, n1.Z);
			Vector3 N2(n2.X, n2.Y, n2.Z);
			Vector3 N3(n3.X, n3.Y, n3.Z);

			float k = -3.0;
			Vector3 N = (N1 + N2 + N3) / k;
			N.normalize();
			Vector3 A(curMesh.Vertices[a].Position.X, curMesh.Vertices[a].Position.Y, curMesh.Vertices[a].Position.Z);
			Vector3 B(curMesh.Vertices[b].Position.X, curMesh.Vertices[b].Position.Y, curMesh.Vertices[b].Position.Z);
			Vector3 C(curMesh.Vertices[c].Position.X, curMesh.Vertices[c].Position.Y, curMesh.Vertices[c].Position.Z);

			HE_Edge* EdgeAB = new HE_Edge;
			HE_Edge* EdgeBC = new HE_Edge;
			HE_Edge* EdgeCA = new HE_Edge;

			EdgeAB->isUsed = false;
			EdgeBC->isUsed = false;
			EdgeCA->isUsed = false;

			EdgeAB->startPoint = A;
			EdgeAB->endPoint = B;
			EdgeBC->startPoint = B;
			EdgeBC->endPoint = C;
			EdgeCA->startPoint = C;
			EdgeCA->endPoint = A;

			EdgeAB->e_pair = nullptr;
			EdgeBC->e_pair = nullptr;
			EdgeCA->e_pair = nullptr;

			EdgeAB->e_succ = nullptr;
			EdgeBC->e_succ = nullptr;
			EdgeCA->e_succ = nullptr;

			EdgeAB->e_succ = EdgeBC;
			EdgeBC->e_succ = EdgeCA;
			EdgeCA->e_succ = EdgeAB;

			if (dualList.empty()) {
				dualList.push_back(EdgeAB);
				dualList.push_back(EdgeBC);
				dualList.push_back(EdgeCA);
			}
			else {
				for (auto iter = dualList.begin(); iter != dualList.end();) {
					if (IsDualEdge(*iter, EdgeAB)) {
						EdgeAB->e_pair = *iter;
						(*iter)->e_pair = EdgeAB;
						iter = dualList.erase(iter);

						continue;
					}

					if (IsDualEdge(*iter, EdgeBC)) {
						EdgeBC->e_pair = *iter;
						(*iter)->e_pair = EdgeBC;
						iter = dualList.erase(iter);
						continue;
					}

					if (IsDualEdge(*iter, EdgeCA)) {
						EdgeCA->e_pair = *iter;
						(*iter)->e_pair = EdgeCA;
						iter = dualList.erase(iter);
						//cout << "CA" << endl;
						continue;
					}
					iter++;
				}

				if (EdgeAB->e_pair == nullptr) {
					dualList.push_back(EdgeAB);
				}
				if (EdgeBC->e_pair == nullptr) {
					dualList.push_back(EdgeBC);
				}
				if (EdgeCA->e_pair == nullptr) {
					dualList.push_back(EdgeCA);
				}

			}

			HE_Face* temp_Face = new HE_Face;
			temp_Face->f_verts.push_back(A);
			temp_Face->f_verts.push_back(B);
			temp_Face->f_verts.push_back(C);

			temp_Face->f_edges.push_back(EdgeAB);
			temp_Face->f_edges.push_back(EdgeBC);
			temp_Face->f_edges.push_back(EdgeCA);

			temp_Face->isGenerated = false;
			temp_Face->isUsed = false;
			temp_Face->isSunkFace = 0;

			temp_Face->normal = N;
			EdgeAB->e_face = temp_Face;
			EdgeBC->e_face = temp_Face;
			EdgeCA->e_face = temp_Face;

			this->face.push_back(temp_Face);
			this->edge.push_back(EdgeAB);
			this->edge.push_back(EdgeBC);
			this->edge.push_back(EdgeCA);

			std::cout.width(7);//i的输出为3位宽
			std::cout << curMesh.Indices.size() - 3 - j;
			std::cout << "\b\b\b\b\b\b\b";

		}
		cout << endl;
		cout << "半边结构生成完毕" << endl;


	}
	return true;
}

bool HE_Mesh::IsFaceUesd() {
	auto e1 = this->edge[0]->isUsed;
	auto e2 = this->edge[1]->isUsed;
	auto e3 = this->edge[2]->isUsed;
	if (e1 == true || e2 == true || e3 == true) {
		return true;
	}

	return false;
}

void HE_Mesh::ResetEdgeState() {
	for (auto iter = this->edge.begin(); iter != this->edge.end(); iter++) {
		(*iter)->isUsed = false;
	}
}

void HE_Mesh::ResetFaceState() {
	for (auto iter = this->face.begin(); iter != this->face.end(); iter++) {
		(*iter)->isGenerated = false;
	}
}

void HE_Mesh::addSunkFace(HE_Face* f) {
	f->isGenerated = false;

	for (auto iter = f->f_edges.begin(); iter != f->f_edges.end(); iter++) {
		(*iter)->e_pair = nullptr;
		(*iter)->isUsed = false;
	}

	sunkFace.push_back(f);
}

void drawFromeMesh(vector<HE_Mesh*>meshList, float*& cavityVe, int& size) {

	float sizes = 0;
	for (auto iter = meshList.begin(); iter != meshList.end(); iter++) {
		auto faces = (*iter)->face;
		for (auto faceIter = faces.begin(); faceIter != faces.end(); faceIter++) {
			++sizes;
		}
	}

	size = sizes;

	auto newArray = new float[sizes * 18];
	cavityVe = newArray;
	int i = 0;

	for (auto iter = meshList.begin(); iter != meshList.end(); iter++) {
		auto faces = (*iter)->face;

		for (auto faceIter = faces.begin(); faceIter != faces.end(); faceIter++) {
			auto normal = (*faceIter)->normal;
			auto vertex = (*faceIter)->f_verts;

			cavityVe[i + 0] = vertex[0].x; cavityVe[i + 1] = vertex[0].y; cavityVe[i + 2] = vertex[0].z; cavityVe[i + 3] = normal.x; cavityVe[i + 4] = normal.y; cavityVe[i + 5] = normal.z;
			cavityVe[i + 6] = vertex[1].x; cavityVe[i + 7] = vertex[1].y; cavityVe[i + 8] = vertex[1].z; cavityVe[i + 9] = normal.x; cavityVe[i + 10] = normal.y; cavityVe[i + 11] = normal.z;
			cavityVe[i + 12] = vertex[2].x; cavityVe[i + 13] = vertex[2].y; cavityVe[i + 14] = vertex[2].z; cavityVe[i + 15] = normal.x; cavityVe[i + 16] = normal.y; cavityVe[i + 17] = normal.z;
			i += 18;

		}
	}
}



void addDataToBuffer(HE_Mesh& mesh, float* a, float*& Ve, int size, float& newSize) {
	auto edge = mesh.edge;
	mesh.ResetEdgeState();
	mesh.ResetFaceState();
	int i = 0;

	for (auto iter = mesh.face.begin(); iter != mesh.face.end(); iter++) {
		auto v = (*iter)->f_verts;
		auto n = (*iter)->normal;

		a[i + 0] = v[0].x; a[i + 1] = v[0].y; a[i + 2] = v[0].z; a[i + 3] = n.x; a[i + 4] = n.y; a[i + 5] = n.z;
		a[i + 6] = v[1].x; a[i + 7] = v[1].y; a[i + 8] = v[1].z; a[i + 9] = n.x; a[i + 10] = n.y; a[i + 11] = n.z;
		a[i + 12] = v[2].x; a[i + 13] = v[2].y; a[i + 14] = v[2].z; a[i + 15] = n.x; a[i + 16] = n.y; a[i + 17] = n.z;
		i += 18;
	}

	i = 0;

	vector<Vector3>normal;
	vector<Vector3>vertex;

	for (auto iter = edge.begin(); iter != edge.end(); iter++) {
		auto currentEdge = (*iter);

		if (currentEdge->isUsed == false) {
			currentEdge->isUsed = true;
			auto dualEdge = currentEdge->e_pair;
			if (dualEdge != nullptr) {
				dualEdge->isUsed = true;
				auto currentFace = currentEdge->e_face;
				auto dualFace = dualEdge->e_face;
				//判断两个三角形是否是凹陷
				if (IsSunk(currentFace, dualFace, currentEdge) == true) {

					if (currentFace->isGenerated == false) {

						bool state = true;
						if (state == true) {
							currentFace->isGenerated = true;
							currentFace->isSunkFace = 1;

							auto v = currentFace->f_verts;
							auto n = currentFace->normal;

							vertex.push_back(v[0]);
							vertex.push_back(v[1]);
							vertex.push_back(v[2]);

							normal.push_back(n);
						}

					}
					if (dualFace->isGenerated == false) {

						bool state = true;

						if (state == true) {
							dualFace->isGenerated = true;
							dualFace->isSunkFace = 1;

							auto v = dualFace->f_verts;
							auto n = dualFace->normal;

							vertex.push_back(v[0]);
							vertex.push_back(v[1]);
							vertex.push_back(v[2]);

							normal.push_back(n);
						}



					}
				}

			}
		}
	}


	newSize = normal.size();
	auto p = new float[newSize * 18];
	Ve = p;
	i = 0;
	int j = 0;
	int k = 0;

	for (auto iter = normal.begin(); iter != normal.end(); iter++) {

		Ve[i + 0] = vertex[k].x; Ve[i + 1] = vertex[k].y; Ve[i + 2] = vertex[k].z; Ve[i + 3] = normal[j].x; Ve[i + 4] = normal[j].y; Ve[i + 5] = normal[j].z;
		Ve[i + 6] = vertex[k + 1].x; Ve[i + 7] = vertex[k + 1].y; Ve[i + 8] = vertex[k + 1].z; Ve[i + 9] = normal[j].x; Ve[i + 10] = normal[j].y; Ve[i + 11] = normal[j].z;
		Ve[i + 12] = vertex[k + 2].x; Ve[i + 13] = vertex[k + 2].y; Ve[i + 14] = vertex[k + 2].z; Ve[i + 15] = normal[j].x; Ve[i + 16] = normal[j].y; Ve[i + 17] = normal[j].z;
		i += 18;
		j += 1;
		k += 3;
	}
	//isSunkFace = 1为内部凹陷，外部凸起，这里需要的是外部凹陷，即内部凸起
	for (auto iter = mesh.face.begin(); iter != mesh.face.end(); iter++) {
		if ((*iter)->isSunkFace == 0) {
			mesh.addSunkFace(*iter);
		}
	}


}

HE_Mesh* createHeMesh(vector<HE_Face*>& face) {
	list<HE_Edge*>dualList;
	HE_Mesh* mesh = new HE_Mesh;
	for (auto faceIter = face.begin(); faceIter != face.end(); faceIter++) {
		auto N = (*faceIter)->normal;
		//法线取反

		//(*faceIter)->normal = (*faceIter)->normal * -1.0;
		auto EdgeAB = (*faceIter)->f_edges[0];
		auto EdgeBC = (*faceIter)->f_edges[1];
		auto EdgeCA = (*faceIter)->f_edges[2];

		if (dualList.empty()) {
			dualList.push_back(EdgeAB);
			dualList.push_back(EdgeBC);
			dualList.push_back(EdgeCA);
		}
		else {
			for (auto dualIter = dualList.begin(); dualIter != dualList.end();) {
				if (IsDualEdge(*dualIter, EdgeAB)) {
					EdgeAB->e_pair = *dualIter;
					(*dualIter)->e_pair = EdgeAB;
					dualIter = dualList.erase(dualIter);

					continue;
				}

				if (IsDualEdge(*dualIter, EdgeBC)) {
					EdgeBC->e_pair = *dualIter;
					(*dualIter)->e_pair = EdgeBC;
					dualIter = dualList.erase(dualIter);

					continue;
				}

				if (IsDualEdge(*dualIter, EdgeCA)) {
					EdgeCA->e_pair = *dualIter;
					(*dualIter)->e_pair = EdgeCA;
					dualIter = dualList.erase(dualIter);

					continue;
				}

				dualIter++;

			}

			if (EdgeAB->e_pair == nullptr) {
				dualList.push_back(EdgeAB);
			}
			if (EdgeBC->e_pair == nullptr) {
				dualList.push_back(EdgeBC);
			}
			if (EdgeCA->e_pair == nullptr) {
				dualList.push_back(EdgeCA);
			}

		}

		mesh->face.push_back(*faceIter);
		mesh->edge.push_back((*faceIter)->f_edges[0]);
		mesh->edge.push_back((*faceIter)->f_edges[1]);
		mesh->edge.push_back((*faceIter)->f_edges[2]);

	}
	return mesh;
}


void splitCavity(HE_Mesh* mesh, vector<HE_Mesh*>& meshList) {
	auto face = mesh->sunkFace;
	auto face1 = createHeMesh(face);
	auto faces = face1->face;

	for (auto faceIter = faces.begin(); faceIter != faces.end(); faceIter++) {
		//这里对凹陷结构集合进行广度优先遍历

		if ((*faceIter)->isUsed == true) {
			continue;
		}


		stack<HE_Face*> faceStack;
		HE_Mesh* newMesh = new HE_Mesh;
		vector<HE_Face*> tempFaceList;

		faceStack.push(*faceIter);
		(*faceIter)->isUsed = true;

		auto epaireA = (*faceIter)->f_edges[0]->e_pair;
		auto epaireB = (*faceIter)->f_edges[1]->e_pair;
		auto epaireC = (*faceIter)->f_edges[2]->e_pair;

		if (epaireA != nullptr) {

			auto faceA = epaireA->e_face;
			if (faceA->isUsed == false) {
				(*faceIter)->f_edges[0]->e_pair->e_face->isUsed = true;
				faceStack.push(faceA);
			}
		}

		if (epaireB != nullptr) {
			auto faceB = epaireB->e_face;
			if (faceB->isUsed == false) {
				(*faceIter)->f_edges[1]->e_pair->e_face->isUsed = true;
				faceStack.push(faceB);
			}
		}

		if (epaireC != nullptr) {
			auto faceC = epaireC->e_face;
			if (faceC->isUsed == false) {
				(*faceIter)->f_edges[2]->e_pair->e_face->isUsed = true;
				faceStack.push(faceC);
			}
		}

		while (!faceStack.empty()) {

			auto currentFace = faceStack.top();
			faceStack.pop();

			auto epaireAs = currentFace->f_edges[0]->e_pair;
			auto epaireBs = currentFace->f_edges[1]->e_pair;
			auto epaireCs = currentFace->f_edges[2]->e_pair;

			if (epaireAs != nullptr) {

				auto faceAs = epaireAs->e_face;
				if (faceAs->isUsed == false) {
					currentFace->f_edges[0]->e_pair->e_face->isUsed = true;
					faceStack.push(faceAs);
				}
			}

			if (epaireBs != nullptr) {

				auto faceBs = epaireBs->e_face;
				if (faceBs->isUsed == false) {
					currentFace->f_edges[1]->e_pair->e_face->isUsed = true;
					faceStack.push(faceBs);
				}
			}

			if (epaireCs != nullptr) {

				auto faceCs = epaireCs->e_face;
				if (faceCs->isUsed == false) {
					currentFace->f_edges[2]->e_pair->e_face->isUsed = true;
					faceStack.push(faceCs);
				}
			}

			tempFaceList.push_back(currentFace);
		}


		//生成半边结构
		auto temp = createHeMesh(tempFaceList);
		meshList.push_back(temp);
	}

}
//对提取出来的空腔进行第一次过滤，按照三角形数量，小于minFaceNum的将被过滤掉
void firstFilteringCavity(vector<HE_Mesh*>& meshList, vector<HE_Mesh*>& newMeshList) {

	int minFaceNum = 10;
	for (auto iter = meshList.begin(); iter != meshList.end(); iter++) {
		if ((*iter)->face.size() > minFaceNum) {
			newMeshList.push_back(*iter);
		}
	}
}

//获取空腔位点位置，具体做法：对于某个空腔，任取其中一个三角形，求三角形中点
//把法线所在的射线移到中点上
void getCavityPosition(vector<HE_Mesh*>& meshList, vector<Vector3>& positions) {
	float k = 1;
	for (auto iter = meshList.begin(); iter != meshList.end(); iter++) {
		auto currentFace = (*iter)->face[0];
		auto normal = currentFace->normal;
		auto midPoint = getMidPoint(currentFace);
		auto result = normal * -k + midPoint;
		positions.push_back(result);
	}
}



vector<string> split(const string& str, const string& pattern)
{
	vector<string> ret;
	if (pattern.empty()) return ret;
	size_t start = 0, index = str.find_first_of(pattern, 0);
	while (index != str.npos)
	{
		if (start != index)
			ret.push_back(str.substr(start, index - start));
		start = index + 1;
		index = str.find_first_of(pattern, start);
	}
	if (!str.substr(start).empty())
		ret.push_back(str.substr(start));
	return ret;
}


void calculateDeepCenter() {

}
