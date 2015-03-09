#include "StdAfx.h"
#include "JointCRF.h"
#include "DenseCRF.h"

string JointCRF::objects[nCLASS] = {"void", "wall", "floor", "picture", "cabinet", // 0-4
	"chair", "table", "window", "door", "ceiling",
	"box", "lamp", "bag", "counter", "bed",  // 10-14
	"books", "blinds", "shelves", "clothes", "bookshelf",
	"desk", "sink", "faucet", "bin", "curtain", // 20-24
	"towel", "mirror", "bowl", "television", "monitor",
	"plant", "jar", "refrigerator", "stove", "printer", // 30-34
	"toilet", "container", "candlestick", "oven", "bathtub",
	"keyboard", "computer" // 40, 41
};

string JointCRF::attributes[nATTRI] = {"wood", "painted", "cotton", "paper", "glass", // 0-4
	"brick", "metal", "leather", "wire", "railing", 
	"plastic", "dirty", "moist-damp", "glossy", "shiny" // 10-14
};


void JointCRF::autoInference(CStr &inDir, CStr &nameNE, CStr &outDir, int borderCrop)
{
	// prepare for global information and load image
	CStr objUnaryDir = inDir + "objUnary/";
	CStr imgDir = inDir + "Images/";
	CStr attUnaryDir = inDir + "AttributeResults/";
	CStr outDataDir = outDir + "Data/";
	CmFile::MkDir(outDataDir);
	Mat img3u = imread(inDir + "Images/" + nameNE + ".jpg"), imgRGB3u;
	CV_Assert(img3u.data != NULL && img3u.isContinuous());
	cvtColor(img3u, imgRGB3u, CV_BGR2RGB);
	int h = img3u.rows, w = img3u.cols;
	Rect reg(borderCrop, borderCrop, w - 2*borderCrop, h - 2*borderCrop);
	CmFile::Copy(inDir + "userinputmasks/" + nameNE + "_P0.png", outDataDir + nameNE + "_UM.png");

	// Load unary potentials, data cost, pairwise correlations
	Mat objUnary1nf = readFileALE(objUnaryDir + nameNE + ".dns");
	vector<Mat> attUnar2f(nATTRI);
	for (int i = 0; i < nATTRI; i++)
		attUnar2f[i] = readFileALE(attUnaryDir + format("/%d/DenseXX/%s.dns", i+1, _S(nameNE)));
	Mat attUnary2nf, objCost1nf, attCost2nf;
	merge(attUnar2f, attUnary2nf);	
	unary2Cost(objUnary1nf, objCost1nf);
	unary2Cost(attUnary2nf, attCost2nf, nATTRI);
	Mat pairAtt1f = readFileKyle(inDir + "att_pairwise.TXT", nATTRI, nATTRI);
	Mat pairAttObj1f = readFileKyle(inDir + "att_joint_pairwise.TXT", nATTRI, nCLASS);

	// Objects and attributes joint segmentation
	CmTimer allTime("AttObjSeg auto");
	allTime.Start();
	Mat objLabel3u, objLablel1u;
	vecM attLabel1u;
	segmentObjAtt(imgRGB3u, objCost1nf, attCost2nf, pairAtt1f, pairAttObj1f, 1, 1, objLablel1u, attLabel1u);
	allTime.StopAndReport();

	// Save automatic inferring results
	imwrite(outDir + nameNE + ".jpg", img3u(reg));
	ActionData::_colorIdx.idx2Color(objLablel1u, objLabel3u);
	imwrite(outDataDir + nameNE + "_OBJ.png", objLablel1u(reg));
	imwrite(outDataDir + nameNE + "_OBJA.png", objLabel3u(reg));
	vecI atts, objs, objsCount(nCLASS);
	vecS attNames, objNames;
	for (int i = 0; i < min(11, nATTRI); i++){
		if (i > 7 && i < 10)
			continue; // Only consider material attributes by now
		if (sum(attLabel1u[i](reg)).val[0] > 255*100)
			atts.push_back(i), attNames.push_back(attributes[i]);
	}
	CV_Assert(objLablel1u.isContinuous());
	byte* label = objLablel1u.data;
	for (int i = 0; i < w*h; i++)
		objsCount[label[i]]++;
	for (int i = 0; i < nCLASS; i++)
		if (objsCount[i] > 100)
			objs.push_back(i), objNames.push_back(objects[i]);
	FileStorage fs(outDataDir + nameNE + "_Inf.yml", FileStorage::WRITE);
	fs<<"ObjNames"<<objNames<<"AttNames"<<attNames;
	Mat objUnarySelect1nf, attUnarySelect2nf, pairAttSel1f, pairAttObjSel1f;
	selectChannel(CmCv::getContinouse(objUnary1nf(reg)), objUnarySelect1nf, objs);
	selectChannel(CmCv::getContinouse(attUnary2nf(reg)), attUnarySelect2nf, atts, 2);
	selectPos(CmCv::getContinouse(pairAtt1f), pairAttSel1f, atts, atts);
	selectPos(CmCv::getContinouse(pairAttObj1f), pairAttObjSel1f, atts, objs);
	CV_Assert(CmCv::writeMat(outDataDir + nameNE + "_ObjS.mat", objUnarySelect1nf));
	CV_Assert(CmCv::writeMat(outDataDir + nameNE + "_AttS.mat", attUnarySelect2nf));
	CV_Assert(CmCv::writeMat(outDataDir + nameNE + "_AttPairS.mat", pairAttSel1f));
	CV_Assert(CmCv::writeMat(outDataDir + nameNE + "_AttObjPairS.mat", pairAttObjSel1f));
}

void JointCRF::interInference(CStr &dir, CStr &nameNE, CMat &sal1f, const ObjDes &des)
{
	CStr dataDir = dir + "Data/";
	FileStorage fs(dataDir + nameNE + "_Inf.yml", FileStorage::READ);
	vecS objNames, attNames;
	fs["ObjNames"] >> objNames;
	fs["AttNames"] >> attNames;
	vecI objs(objNames.size()), atts(attNames.size());
	for(size_t i = 0; i < objs.size(); i++)
		objs[i] = ActionData::findFromList(objNames[i], ActionData::_objects);
	for(size_t i = 0; i < atts.size(); i++)
		atts[i] = ActionData::findFromList(attNames[i], ActionData::_materials);

	Mat objUnary1nf, attUnary2nf, pairAtt1f, pairAttObj1f;
	CV_Assert(CmCv::readMat(dataDir + nameNE + "_OBJS.mat", objUnary1nf));
	CV_Assert(CmCv::readMat(dataDir + nameNE + "_AttS.mat", attUnary2nf));
	CV_Assert(CmCv::readMat(dataDir + nameNE + "_AttPairS.mat", pairAtt1f));
	CV_Assert(CmCv::readMat(dataDir + nameNE + "_AttObjPairS.mat", pairAttObj1f));
	Mat img3u = imread(dir + nameNE + ".jpg"), imgRGB3u;
	cvtColor(img3u, imgRGB3u, CV_BGR2RGB);
	CV_Assert(atts.size()*2 == attUnary2nf.channels());

	// Interaction
	int crntObj = ActionData::findFromList<int>(des.obj, objs);
	if (crntObj != -1)
		add2Chan(objUnary1nf, sal1f, crntObj, 5);
	int crntAtt = ActionData::findFromList(des.mat, atts);
	if (crntAtt != -1){
		add2Chan(attUnary2nf, sal1f, crntAtt*2, 5);
		add2Chan(attUnary2nf, sal1f, crntAtt*2+1, -5);
	}
	
	// Objects and attributes joint segmentation
	Mat objCost1nf, attCost2nf;
	unary2Cost(objUnary1nf, objCost1nf);
	unary2Cost(attUnary2nf, attCost2nf, atts.size());
	Mat objLabel3u, objLablel1u;
	vecM attLabel1u;
	segmentObjAtt(imgRGB3u, objCost1nf, attCost2nf, pairAtt1f, pairAttObj1f, 1, 1, objLablel1u, attLabel1u);

	// Save automatic inferring results
	reLabel(objLablel1u, objs);
	ActionData::_colorIdx.idx2Color(objLablel1u, objLabel3u);
	imwrite(dataDir + nameNE + "_OBJC.png", objLabel3u);
	imwrite(dataDir + nameNE + "_OBJ.png", objLablel1u);
	for (size_t i = 0; i < attLabel1u.size(); i++){
		if (sum(attLabel1u[i]).val[0] > 255*100)
			imwrite(dataDir + nameNE + format("_ATT_%d.png", atts[i]), attLabel1u[i]);
	}
}

void JointCRF::add2Chan(Mat &tagNf, CMat &val1f, int tagC, float scale)
{
	CV_Assert(tagNf.size == val1f.size && tagNf.isContinuous() && val1f.isContinuous());
	const int pixNum = tagNf.cols * tagNf.rows, nChal = tagNf.channels();
	const float *val = (float*)(val1f.data);
	float *tag = (float*)(tagNf.data) + tagC;
#pragma omp parallel for
	for (int p = 0; p < pixNum; p++)
		tag[p*nChal] += scale * val[p];
}

void JointCRF::reLabel(Mat &label1u, const vecI& newLabel)
{
	int rows = label1u.rows, cols = label1u.cols;
	if (label1u.isContinuous())
		rows = 1, cols *= label1u.rows;
	for (int r = 0; r < rows; r++){
		byte* label = label1u.ptr<byte>(r);
#pragma omp parallel for
		for (int c = 0; c < cols; c++)
			label[c] = newLabel[label[c]];
	}
}

void JointCRF::selectChannel(CMat &unariesNf, Mat &sUnariesNf, const vecI &channs, int numPerChan)
{
	CV_Assert(unariesNf.isContinuous() && (unariesNf.channels() % numPerChan == 0));
	const int numPix = unariesNf.rows * unariesNf.cols, slectDim = channs.size() * numPerChan;
	vecI pos;
	for (size_t i = 0; i < channs.size(); i++){
		CV_Assert(channs[i] * numPerChan < unariesNf.channels());
		for (int c = 0; c < numPerChan; c++)
			pos.push_back(channs[i]*numPerChan + c);
	}
	sUnariesNf = Mat(unariesNf.size(), CV_MAKETYPE(CV_32F, slectDim));
#pragma omp parallel for
	for (int i = 0; i < numPix; i++){
		const float* src = (float*)(unariesNf.data) + i * unariesNf.channels();
		float *dst = (float*)(sUnariesNf.data) + i * sUnariesNf.channels();
		for (int d = 0; d < slectDim; d++)
			dst[d] = src[pos[d]];
	}
}

void JointCRF::selectPos(CMat &mat1f, Mat &matSel1f, const vecI &rows, const vecI &cols)
{
	matSel1f = Mat(rows.size(), cols.size(), CV_32F);
	for (size_t r = 0; r < rows.size(); r++){
		const float* src = mat1f.ptr<float>(rows[r]);
		float* dst = matSel1f.ptr<float>(r);
		for (size_t c = 0; c < cols.size(); c++)
			dst[c] = src[cols[c]];
	}
}

void JointCRF::segmentObjAtt(CMat &img3u, CMat &objCost1nf, CMat &attCost2nf, CMat &pairAtt, 
	CMat &pairAttObj, float wJoin, float wCor, Mat &objLablel1s, vecM &attLabel1s)
{
	int nAtt = pairAttObj.rows, nObj = pairAttObj.cols;
	DenseCRF2D crfPlane(img3u.cols, img3u.rows, nObj, nAtt);
	crfPlane.setUnaryEnergy((float*)objCost1nf.data);
	crfPlane.addPairwiseGaussian(3, 3, 3);
	crfPlane.addPairwiseBilateral(50, 50, 15, 15, 15, img3u.data, 5);
	attLabel1s.resize(nAtt);

	CmTimer allTime("AttObjSeg");
	allTime.Start();
	crfPlane.set_attributes((float*)attCost2nf.data, (float*)pairAtt.data, (float*)pairAttObj.data, wJoin, wCor);
	objLablel1s = crfPlane.map(5, attLabel1s);
	allTime.StopAndReport();
}

// each unaryNf.channels()/nGroup channels are normalized to get the cost
void JointCRF::unary2Cost(CMat &unaryNf, Mat &costNf, int nGroup)
{
	const int tChan = unaryNf.channels(); // Total number of channels
	const int nChan = tChan/nGroup; // Number of channels within which needs to be normalized
	const int pixNum = unaryNf.cols * unaryNf.rows;
	unaryNf.copyTo(costNf);
	CV_Assert(costNf.isContinuous() && (tChan % nGroup == 0));
#pragma omp parallel for
	for (int p = 0; p < pixNum; p++){
		float *res = (float*)(costNf.data) + p*tChan;
		for (int i = 0; i < nGroup; i++, res += nChan){
			float sum = 1e-8f;
			for (int c = 0; c < nChan; c++)
				sum += (res[c] = exp(res[c]));
			for (int c = 0; c < nChan; c++)
				res[c] = -log(res[c]/sum);
		}
	}
}

// return float type and constrain data range to [-10, 10]
Mat JointCRF::readFileALE(CStr &dnsFileName)
{
	FILE *fp = fopen(_S(dnsFileName), "rb");
	CV_Assert_(fp != NULL, ("Can't load file %s\n", _S(dnsFileName)));
	int temp[3]; 
	fread(temp, sizeof(int), 3, fp);
	Mat dataM(temp[1], temp[0], CV_MAKETYPE(CV_32F, temp[2]));
	int dataNum = temp[0]*temp[1]*temp[2];
	double *dataBuf = new double[dataNum];
	float *fData = (float*)dataM.data;
	fread(dataBuf, sizeof(double), dataNum, fp);
#pragma omp parallel for
	for (int i = 0; i < dataNum; i++)
		fData[i] = max(-10.0, min(dataBuf[i], 10.0));
	delete []dataBuf;
	fclose(fp);
	flip(dataM, dataM, CV_FLIP_VERTICAL);
	return dataM;
}

// Read file's from Kyle
Mat JointCRF::readFileKyle(CStr &fileName, int height, int width)
{
	FILE *f = fopen(_S(fileName), "rb");
	Mat mat(height, width, CV_32F);
	fread(mat.data, sizeof(float), height * width, f);
	fclose(f);
	return mat;
}

// D:\jointobjatt3fast\ 68.jpg D:\WkDir\ImgSpirit\ 10 4 0
int JointCRF::demo(int argc, char* argv[])
{
	if (argc != 7){
		printf("Input command error at %d:%s\n", __LINE__, __FILE__);
		return 1;
	}

	CStr dir = argv[3], nameNE = CmFile::GetNameNE(argv[2]), inDir = argv[1];
	ActionData::loadGlobalSettings("Setting.yml");
	if (!CmFile::FileExist(dir + "Data/" + nameNE + "_OBJS.mat"))
		autoInference(argv[1], nameNE, dir, atoi(argv[4]));
	Mat sal1f = imread(dir + "Data/" + nameNE + "_UM.png", CV_LOAD_IMAGE_GRAYSCALE);
	sal1f.convertTo(sal1f, CV_32F, 2.0/255, -1);
	CV_Assert(sal1f.data != NULL);
	ObjDes des;
	des.obj = atoi(argv[5]), des.mat = atoi(argv[6]);
	interInference(dir, nameNE, sal1f, des);
	return 0;
}