#include <SOP/SOP_Node.h>
#include <GEO/GEO_Hedge.h>
#include <GEO/GEO_HedgeInterface.h>
#include <GEO/GEO_PolyInterface.h>
#include <OP/OP_AutoLockInputs.h>
#include <GA/GA_Handle.h>
#include <GA/GA_Types.h>
#include <OP/OP_OperatorTable.h>
#include <OP/OP_SaveFlags.h>
#include <PRM/PRM_Include.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class GEO_PolyInterface;

class DIFF_SOP : public SOP_Node {


public:
	typedef  Eigen::Matrix<fpreal, Eigen::Dynamic, Eigen::Dynamic> EigenMat;
	typedef  Eigen::SparseMatrix<fpreal > SparseMat;
	typedef  Eigen::Triplet<fpreal> coeffs;

	static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);
	static PRM_Template  MyTemplateList[];

	UT_Array<UT_Vector2> edgeList;
	GEO_PolyInterface *polinterface;
	GEO_HedgeInterface *hedgeinterface;

private:

	fpreal cotangent(GEO_Hedge hedge);
	void gradient(GA_ROHandleF scalarfield, GA_Offset primOffset);
	GA_ROHandleF  sclar;
	void Laplacian(GA_ROHandleV3 vectorfield, OP_Context &context);
	void divergence(GA_ROHandleV3 vectorfield, GA_Offset primOffset);
	GA_ROHandleV3  vectorfield;
    void Vectordecomposition(GA_ROHandleV3 vectorfield, OP_Context &context);
	UT_Vector3 oneformtofield(UT_Array<fpreal> &oneformlist, GA_Offset primoffset);
	void covector(UT_Array<fpreal> &oneformlist, GA_ROHandleV3 &vectorfield);

	SparseMat deriative0(SparseMat d0);
	SparseMat deriative1(SparseMat d1);
	SparseMat star1(SparseMat star1);
	SparseMat star2(SparseMat star);
	SparseMat star2inv(SparseMat star);
	SparseMat star3(SparseMat);

	int getedgeix(GEO_Hedge h);

	UT_Array<int>  primhedgelist;
	UT_Array<int>  primhedgelistbnd;

	int    fieldchoice();
	void SCALAR(UT_String &str, fpreal t)
	{
		evalString(str, "Field", 0, t);
	}


protected:
	DIFF_SOP(OP_Network * net, const char *name, OP_Operator *op);
	virtual ~DIFF_SOP();
	virtual OP_ERROR cookMySop(OP_Context &context);





};


