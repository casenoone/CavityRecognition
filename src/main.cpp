#include <iostream>
#include "HalfEdge.h"
#include "vector3.h"
#include <iomanip>
#include <sstream>
using namespace std;


//��������ļ�����0000.obj  0001.obj ...����������main.cppͬ���ļ���
//�޸�n�Ĵ�С���������һ���ļ���0095,���ڴ�0000��ʼ��������ôn��=95

//���Ҫ�ٴμ��㣬Ҫ��ɾ���ļ������Ѿ���õ�txt�ļ�

//�ó���ʹ���˵������⣺obj_loader
//����������û��дcpp�ļ�������HalfEdge.h�����obj_loader.h�����
//��������ҵ����ض���....��
//�ڽ���������������Ҽ���Ŀ����-������-����������-��������-������ѡ��-������ѡ��
//���������  /force
int main()
{

    objl::Loader loader;

    string name;
    int n = 96;
    for (int i = 0; i < n; ++i) {
        
        name = "";
        string temp = to_string(i);
        int tempK = i / 10;

        //һλ��
        if (tempK == 0) {
            name = "000" + temp;
        }

        //��λ��
        else if (tempK >= 1&&tempK<10) {
            name = "00" + temp;
        }

        //��λ��
        else if(tempK >= 10 && tempK < 100) {
            name = "0" + temp;
        }
        //��λ��
        else {

        }
        
        name = name;
        
        cout << name << endl;

        HE_Mesh* mesh = new HE_Mesh;
        mesh->LoadFromObj(loader, name+".obj");
        int size = mesh->face.size();

        float* preVe = new float[size * 18];
        float* Ve = nullptr;

        float sizeofs = size * 18 * 4;
        float newSizeof;
        addDataToBuffer(*mesh, preVe, Ve, size, newSizeof);

        auto l1 = mesh->sunkFace.size();

        vector<Vector3> cavityPositions;

        vector<HE_Mesh*> tempMeshList;
        vector<HE_Mesh*> meshList;
        splitCavity(mesh, tempMeshList);
        firstFilteringCavity(tempMeshList, meshList);
        getCavityPosition(meshList, cavityPositions);
        cout << "��ǻ��Ŀ��" << meshList.size() << endl;



        ofstream out(name+".txt" , ios::app);
        for (auto iter = cavityPositions.begin(); iter != cavityPositions.end(); iter++) {

            auto x = iter->x;
            auto y = iter->y;
            auto z = iter->z;
            out << x << "," << y << "," << z << endl;
        }


        delete mesh;
    }

    

    return 0;
}





