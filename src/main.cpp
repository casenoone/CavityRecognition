#include <iostream>
#include "HalfEdge.h"
#include "vector3.h"
#include <iomanip>
#include <sstream>
using namespace std;


//将待算的文件按照0000.obj  0001.obj ...命名，放入main.cpp同级文件夹
//修改n的大小，例如最后一个文件是0095,由于从0000开始命名，那么n就=95

//如果要再次计算，要先删掉文件夹里已经算好的txt文件

//该程序使用了第三方库：obj_loader
//可能由于其没有写cpp文件，于是HalfEdge.h里包含obj_loader.h会出错
//如果报错“找到多重定义....”
//在解决方案管理器里右键项目名称-》属性-》配置属性-》链接器-》所有选项-》附加选项
//在里面添加  /force
int main()
{

    objl::Loader loader;

    string name;
    int n = 96;
    for (int i = 0; i < n; ++i) {
        
        name = "";
        string temp = to_string(i);
        int tempK = i / 10;

        //一位数
        if (tempK == 0) {
            name = "000" + temp;
        }

        //两位数
        else if (tempK >= 1&&tempK<10) {
            name = "00" + temp;
        }

        //三位数
        else if(tempK >= 10 && tempK < 100) {
            name = "0" + temp;
        }
        //四位数
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
        cout << "空腔数目：" << meshList.size() << endl;



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





