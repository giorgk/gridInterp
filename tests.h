#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "gridInterp.h"

void test1D(){

    // Read the data for the POINT mode and LINEAR method
    GRID_INTERP::interp<1> oneD;
    oneD.getDataFromFile("Rgridinterp/test1_data_p.tmp");
    
    std::vector<double> vv;
    std::vector<std::vector<double> > xx;
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1D_pl.tmp", xx, vv);

    // For the POINT mode and NEAREST method we use the same data but set different methods

    oneD.setModeMethod(GRID_INTERP::MODE::POINT, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1D_pn.tmp", xx, vv);


    // For the CELL mode the number of axis data are one more than the values
    oneD.reset();
    oneD.getDataFromFile("Rgridinterp/test1_data_c.tmp");

    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1D_cl.tmp", xx, vv);

    oneD.setModeMethod(GRID_INTERP::MODE::POINT, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1D_cn.tmp", xx, vv);

}
