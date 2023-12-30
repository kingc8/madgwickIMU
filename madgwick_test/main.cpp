//
//  main.cpp
//  madgwickFilter test.
//
//  Created by Blake Johnson on 4/28/20.
//  Additions by Cameron King on 29/12/23.
//

#include <iostream>
#include "madgwickIMU.hpp"

using namespace std;

int main()
{
    madgwick mw;

    quat q = mw.MadgwickAHRSupdateIMU({ 0.05, 0.065, 0.9 }, { 0.0, 0.0, 0.0 });
    cout << "Quarternion = " << q.x << ", " << q.y << ", " << q.z <<", " << q.s << "\n";

    return 0;
}
