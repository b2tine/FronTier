#ifndef DRAG_H
#define DRAG_H
#include <algorithm>
#include "spring_solver.h"

class SpringVertex;

class Drag {
    void operator=(const Drag&);
    Drag(const Drag&);
    static double m_thickness;
    static double m_tol;
public:
    double m_t;
    
    //set point type and register points
    virtual void preprocess(std::vector<SpringVertex*>&){}; 

    //update parameters for next step
    virtual void postprocess(std::vector<SpringVertex*>&){}; 

    //set velocity for each spring vertex
    //according to point type
    virtual void setVel(SpringVertex*) = 0;
    virtual void updateVel(std::vector<SpringVertex*>& pts, double t) {}

    //set acceleration for each spring vertex
    //according to point type
    virtual void setAccel(SpringVertex*) = 0;

    //return number of parameters needed
    virtual int getParameterNumber() {return 0;}

    void setTimeStepSize(double dt) {m_dt = dt;}
    double getTimeStepSize() {return m_dt;}
    void spinToAxis(double*, double*, double, double*); 
    
    static void setThickness(double t){m_thickness = t;}
    static double getThickness() {return m_thickness;}
    static void setTolerance(double t){m_tol = t;}
    static double getTolerance() {return m_tol;}

    //prototype pattern
    //the drag will be created using 
    //dragFactory() by providing appropriate 
    //data, should not call constructor directly
    friend class DragProtoInit;
    class Info {
	std::string _id;
	std::vector<double> _data;
     public:
	Info(std::string ident, std::vector<double> &dat):
		_id(ident), _data(dat) {}
	Info() {}
	std::string &id() {return _id;}
	const std::string& id() const {return _id;}
	std::vector<double>& data() {return _data;}
	const std::vector<double>& data() const {return _data;}
	void clear() {_id.clear(); _data.clear(); }
        bool empty() { return _id.empty(); }
    };
    virtual std::string id() = 0; 
    static std::vector<Drag*> prototypes;
    static Drag* dragFactory(const Info &);
    virtual Drag* clone(const Info &) = 0;
    virtual size_t dataSize() { return 0; }
protected:
    bool validateData(const Drag::Info&); 
    enum {FREE_POINT,ROTATE_POINT, STATIC_POINT, TRANS_DPOINT, 
		TRANS_CPOINT1, TRANS_CPOINT2};
    double m_dt;
    bool first;
    Drag(double t) : m_t(t) { first = true; }
    Drag(){ first = true; }
};

class PLANE
{
public: 
    double center[3];
    double dir[3];
    double nor[3];

    bool isFront(double* p);    
    bool isBack(double* p);
    double distance(double* p);
    PLANE(const double*, const double*);
};

class LINE
{
public:
    double pOnLine1[3]; 
    double pOnLine2[3];

    LINE(const double*, const double*);
    double distance(double* p); 
};

class PointDrag : public Drag {
public:
    double cent[3];
    double rad;
    double m_v[3];
    double m_a[3];
    virtual void preprocess(std::vector<SpringVertex*>&);
    virtual void setVel(SpringVertex*);
    virtual void setAccel(SpringVertex*);

    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "PointDrag";}
    virtual size_t dataSize() {return 11;}
    PointDrag(const double c[], double r, const double v[], const double a[],  double t);
    PointDrag(){};
};

class GravityDrag : public PointDrag {
public:
    void setAccel(SpringVertex*); 
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "GravityDrag";}
    virtual size_t dataSize() {return 8;}
    GravityDrag(const double c[], const double r, const double a[], double t);
    GravityDrag(){}
};

class LineDrag : public Drag {
protected :
    LINE* dragLine;
    LINE* controlLine1; 
    LINE* controlLine2; 
    double veld[3];
    double velc1[3];
    double accelc1[3];
    double velc2[3];
    double accelc2[3];
    double accelStartTime1; // time for starting accelating motion 
    double accelStartTime2; // time for starting accelating motion
    double m_ct;  // current total time for this folding plan
public :
    LineDrag() { m_ct = 0.0; }
    virtual void preprocess(std::vector<SpringVertex*>&); 
    virtual void postprocess(std::vector<SpringVertex*>&);
    virtual void setVel(SpringVertex*); 
    virtual void setAccel(SpringVertex*);  
    virtual Drag* clone(const Drag::Info&); 
    virtual size_t dataSize() {return 36;}
    virtual std::string id() { return "LineDrag"; }
    void accumCurTime() { m_ct += m_dt; }
};

class GravityBoxDrag : public Drag 
{
public:
    double L[3];
    double U[3];
    double g[3];
    void preprocess(std::vector<SpringVertex*>&);
    void setVel(SpringVertex*) {};
    void setAccel(SpringVertex*);
    Drag* clone(const Drag::Info&);
    std::string id() {return "GravityBoxDrag";}
    virtual size_t dataSize() {return 10;}
    GravityBoxDrag(double t);
    GravityBoxDrag(){}
};

class MultiplePointDrag : public Drag {
public:
    struct DragPoint {
	double cent[3];
	double rad;
	double m_v[3];
	double m_a[3];
	DragPoint(double[],double,double[],double[]);
    };
    std::vector<DragPoint*> dragPoints;
    void preprocess(std::vector<SpringVertex*>&);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    int getDragPointNumber(double []);
    void addDragPoint(DragPoint*);
    virtual Drag* clone(const Drag::Info&);
    std::string id() {return "MultiplePointDrag";}
    MultiplePointDrag(double t);
    MultiplePointDrag(){}
};

class FoldDrag : public Drag {
public:
    PLANE* static_plane;
    double angVel;
    double spinOrig[3];
    double spinDir[3];
    
    void preprocess(std::vector<SpringVertex*>&);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*);
    std::string id() {return "FoldDrag";}
    virtual Drag* clone(const Drag::Info&);
    virtual size_t dataSize() { return 14; }
    FoldDrag() {}
};

class ZFoldDrag: public Drag {
private : 
    PLANE* static_plane1;
    PLANE* static_plane2;  
    double angVel; 
    double spinOrig[3]; 
    double spinDir[3]; 
public :
    virtual void preprocess(std::vector<SpringVertex*>&);
    virtual void setVel(SpringVertex*);
    virtual void setAccel(SpringVertex*); 
    virtual std::string id() { return "ZFoldDrag"; }
    virtual Drag* clone(const Drag::Info&); 
    virtual size_t dataSize() { return 20; }
    ZFoldDrag() {}
}; 

class CloseUmbrellaDrag : public Drag {
    int num_tuck_line;
    double spinOrig[3];
    double spinDir[3];
    double angVel;
public:
    void preprocess(std::vector<SpringVertex*>&);
    std::string id() {return "CloseUmbrellaDrag";}
    Drag* clone(const Drag::Info&);
    void setVel(SpringVertex*);
    void setAccel(SpringVertex*){}
    virtual size_t dataSize() { return 9; }
    CloseUmbrellaDrag() {}
};

class RelaxDrag : public Drag {
public: 
    std::string id() {return "RelaxDrag";}
    Drag* clone (const Drag::Info&);
    void preprocess(std::vector<SpringVertex*>&);
    void setVel(SpringVertex* sv) {}
    void setAccel(SpringVertex* sv) {}
};

class CompressDrag: public Drag {
protected:
    double center[3];
    double accel[3];
    double thickness;
public:
    void preprocess(std::vector<SpringVertex*>&);
    std::string id() {return "CompressDrag";}
    Drag* clone(const Drag::Info&);
    void setVel(SpringVertex* sv) {}
    void setAccel(SpringVertex* sv);
    virtual size_t dataSize() { return 8; }
};

class SeparateDrag: public Drag {
    double spin_center[3];
    double spin_dir[3];
    double nor[3];
    double angVel;
    double radius;
    double old_body_center[3];
public:
    void preprocess(std::vector<SpringVertex*>&);
    void postprocess(std::vector<SpringVertex*>&);
    std::string id() {return "SeparateDrag";}
    Drag* clone(const Drag::Info&);
    void setVel(SpringVertex* sv);
    void setAccel(SpringVertex* sv){}
    virtual size_t dataSize() { return 11; }
};


class RollDrag: public Drag {
    double spin_center[3];
    double spin_dir[3];

    double mov_center[3];
    double mov_dir[3];

    int num_layers;
    double ang_vel;
    double spacing;
public:
    void preprocess(std::vector<SpringVertex*>&);
    std::string id() {return "RollDrag";}
    Drag* clone(const Drag::Info&);
    virtual size_t dataSize() {return 14;}
    void setVel(SpringVertex* sv);
    void setAccel(SpringVertex* sv){}
};

class AlignDrag: public Drag {
    double rotate_center[3];
    double rotate_axis[3];
    double gravity_center[3];
    double dir[3];
    double rotate_angle;
    double R[3][3];
    AlignDrag(const double[], const double[]);
public:
    AlignDrag(){ rotate_angle = 0; };
    void preprocess(std::vector<SpringVertex*>&);
    void postprocess(std::vector<SpringVertex*>&);
    std::string id() {return "AlignDrag";}
    Drag* clone(const Drag::Info&);
    virtual size_t dataSize() {return 6;}
    void setVel(SpringVertex* sv){}
    void setAccel(SpringVertex* sv){}
};

#endif
