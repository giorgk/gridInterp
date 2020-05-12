#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "gridInterp.h"


void test2D_const(){

    double fromX = -21.0;
    double toX = 21.0;
    double dx = 0.2;

    double fromY = -10.0;
    double toY = 9.5;
    double dy = 0.2;

    // Read the data for the POINT mode and LINEAR method
    GRID_INTERP::interp<2> testD;
    testD.getDataFromFile("Rgridinterp/test2D_data_var_p.tmp");
    
    std::vector<double> vv;
    std::vector<std::vector<double> > xx;
    for(double y = fromY; y < toY; y = y+dy){
        for(double x = fromX; x < toX; x = x+dx){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(testD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test2D_pl.tmp", xx, vv);

    testD.setModeMethod(GRID_INTERP::MODE::POINT, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double y = fromY; y < toY; y = y+dy){
        for(double x = fromX; x < toX; x = x+dx){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(testD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test2D_pn.tmp", xx, vv);

    testD.reset();
    testD.getDataFromFile("Rgridinterp/test2D_data_var_c.tmp");
    
    vv.clear();
    xx.clear();
    for(double y = fromY; y < toY; y = y+dy){
        for(double x = fromX; x < toX; x = x+dx){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(testD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test2D_cl.tmp", xx, vv);

    testD.setModeMethod(GRID_INTERP::MODE::CELL, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double y = fromY; y < toY; y = y+dy){
        for(double x = fromX; x < toX; x = x+dx){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(testD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test2D_cn.tmp", xx, vv);

    return;
    /*
    // Read the data for the POINT mode and LINEAR method
    GRID_INTERP::interp<2> twoD;
    twoD.getDataFromFile("Rgridinterp/peak_data_cnst_p.tmp");
    
    std::vector<double> vv;
    std::vector<std::vector<double> > xx;
    for(double y = -29.4; y < 30.3; y = y+0.5047){
        for(double x = -29.5; x < 30.5; x = x+0.4856){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(twoD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_peak_pl.tmp", xx, vv);

    // For the POINT mode and NEAREST method we use the same data but set different methods

    twoD.setModeMethod(GRID_INTERP::MODE::POINT, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double y = -29.4; y < 30.3; y = y+0.5047){
        for(double x = -29.5; x < 30.5; x = x+0.4856){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(twoD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_peak_pn.tmp", xx, vv);
    */
}

void test1D_var(){
    // Read the data for the POINT mode and LINEAR method
    GRID_INTERP::interp<1> oneD;
    oneD.getDataFromFile("Rgridinterp/test1_data_var_p.tmp");
    
    std::vector<double> vv;
    std::vector<std::vector<double> > xx;
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test1D_pl.tmp", xx, vv);

    // For the POINT mode and NEAREST method we use the same data but set different methods

    oneD.setModeMethod(GRID_INTERP::MODE::POINT, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test1D_pn.tmp", xx, vv);


    // For the CELL mode the number of axis data are one more than the values
    oneD.reset();
    oneD.getDataFromFile("Rgridinterp/test1_data_var_c.tmp");

    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test1D_cl.tmp", xx, vv);

    oneD.setModeMethod(GRID_INTERP::MODE::CELL, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test1D_cn.tmp", xx, vv);
}

void test1D_const(){

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

    oneD.setModeMethod(GRID_INTERP::MODE::CELL, GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1D_cn.tmp", xx, vv);

}
