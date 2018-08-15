#include "BV_Node.h"

#ifndef BVH_H
#define BVH_H



class BVH
{
    private:
        
        BV_Node* root{nullptr};

    public:

        const BV_Node* const getRoot() const {return root;}



};




#endif
