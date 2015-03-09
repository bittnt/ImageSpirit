#pragma once

/************************************************************************/
/* Attributes based interaction in a JointCRF framework                 */
/* Data files, object names, attributes names, are in ActionData class	*/
/************************************************************************/
#include "ActionData.h"

class JointCRF
{
public:
	static const int nCLASS = 42, nATTRI = 15;
	static string objects[nCLASS];
	static string attributes[nATTRI];

	static void autoInference(CStr &inDir, CStr &nameNE, CStr &outDir, int borderCrop = 0);

	static void interInference(CStr &dir, CStr &nameNE, CMat &sal1f, const ObjDes &des);

	static void selectChannel(CMat &unariesNf, Mat &sUnariesNf, const vecI &channs, int numPerChan = 1);

	static void selectPos(CMat &mat1f, Mat &matSel1f, const vecI &rows, const vecI &cols);

	static void reLabel(Mat &label1u, const vecI& newLabel);

	static void segmentObjAtt(CMat &img3u, CMat &objCost1nf, CMat &attrCost2nf, CMat &pairAtt, 
		CMat &pairAttObj, float wJoin, float wCor, Mat &objLablel1s, vecM &attLabel1s);

	static int demo(int argc, char* argv[]);

public: 
	// return float type and constrain data range to [-10, 10]
	static Mat readFileALE(CStr &dnsFileName); 

	// read data file from Kyle
	static Mat readFileKyle(CStr &fileName, int height, int width);

	// each unaryNf.channels()/nGroup channels are normalized to get the cost
	static void unary2Cost(CMat &unaryNf, Mat &costNf, int nGroup = 1);

	static void loadObjAtrData(CStr &objDir, CStr &atrDir, Mat &objUnar1nf, Mat &atrUnar2nf);	

	static void add2Chan(Mat &tagNf, CMat &val1f, int tagC, float scale = 1);
};

