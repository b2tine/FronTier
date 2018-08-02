#include "drag.h"
#include "origami.h"
#include <vector>

std::vector<Drag*> Drag::prototypes;

//this class need to be modified
//whenever adding a new type of drag
class DragProtoInit {
    
    PointDrag pd;
    LineDrag ld; 
    GravityDrag gd;
    MultiplePointDrag md;
    FoldDrag sd;
    ZFoldDrag zsd; 
    CloseUmbrellaDrag cud;
    RelaxDrag rd;
    GravityBoxDrag gbd;
    CompressDrag cpd;
    SeparateDrag spd;
    RollDrag rold; 
    OrigamiFold ogmd;
    AlignDrag ag;
    DragProtoInit() {
	Drag::prototypes.push_back(&pd);
	Drag::prototypes.push_back(&ld); 
	Drag::prototypes.push_back(&gd);
	Drag::prototypes.push_back(&md);
	Drag::prototypes.push_back(&sd);
        Drag::prototypes.push_back(&zsd);
	Drag::prototypes.push_back(&cud);
	Drag::prototypes.push_back(&rd);
	Drag::prototypes.push_back(&gbd);
	Drag::prototypes.push_back(&cpd);
	Drag::prototypes.push_back(&spd);
        Drag::prototypes.push_back(&rold);
        Drag::prototypes.push_back(&ogmd);
        Drag::prototypes.push_back(&ag);
    }
    static DragProtoInit singleton;
};

DragProtoInit DragProtoInit::singleton;
