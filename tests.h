#pragma once

#include <iostream>
#include <vector>
#include <string>

#include "gridInterp.h"


void test2D(){

    std::cout << "2D Interpolation using evenly spaced data" << std::endl;
    GRID_INTERP::interp<2> peaks;
    peaks.getDataFromFile("Rgridinterp/peak_data_cnst.tmp");
    
    std::vector<double> vv;
    std::vector<std::vector<double> > xx;
    for(double y = -29.42; y < 30.54; y = y+0.5047){
        for(double x = -29.52; x < 30.64; x = x+0.4856){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(peaks.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_peak_l.tmp", xx, vv);
    peaks.setMethod(GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double y = -29.4; y < 30.3; y = y+0.5047){
        for(double x = -29.5; x < 30.5; x = x+0.4856){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(peaks.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_peak_n.tmp", xx, vv);

    std::cout << "2D Interpolation using variably spaced data" << std::endl;

    double fromX = -21.0;
    double toX = 21.0;
    double dx = 0.2;

    double fromY = -10.0;
    double toY = 9.5;
    double dy = 0.2;

    // Read the data for the POINT mode and LINEAR method
    GRID_INTERP::interp<2> testD;
    testD.getDataFromFile("Rgridinterp/test2D_data_var.tmp");
    
    vv.clear();
    xx.clear();
    for(double y = fromY; y < toY; y = y+dy){
        for(double x = fromX; x < toX; x = x+dx){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(testD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test2D_l.tmp", xx, vv);

    testD.setMethod(GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double y = fromY; y < toY; y = y+dy){
        for(double x = fromX; x < toX; x = x+dx){
            xx.push_back(std::vector<double> {x,y});
            vv.push_back(testD.interpolate(x, y));
        }
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test2D_n.tmp", xx, vv);
}

void test1D(){
    std::cout << "1D Interpolation using evenly spaced data" << std::endl;
    // Read the data for LINEAR method 
    GRID_INTERP::interp<1> oneD;
    oneD.getDataFromFile("Rgridinterp/test1_cnst_data.tmp");
    
    std::vector<double> vv;
    std::vector<std::vector<double> > xx;
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1_l.tmp", xx, vv);

    // Set the method to NEAREST and repeat the interpolation

    oneD.setMethod(GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 39.0; x < 66.5; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_cnst_test1_n.tmp", xx, vv);

    std::cout << "1D Interpolation using variably spaced data" << std::endl;
    oneD.reset();
    oneD.getDataFromFile("Rgridinterp/test1_var_data.tmp");
    vv.clear();
    xx.clear();
    for(double x = 38.0; x < 82.2; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test1_l.tmp", xx, vv);

    oneD.setMethod(GRID_INTERP::METHOD::NEAREST);
    vv.clear();
    xx.clear();
    for(double x = 38.0; x < 82.2; x = x+0.1){
        xx.push_back(std::vector<double> {x});
        vv.push_back(oneD.interpolate(x));
    }
    GRID_INTERP::writeCoordsValues("Rgridinterp/res_var_test1_n.tmp", xx, vv);
}

