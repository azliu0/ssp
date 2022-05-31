# testing script for SSP 2021 Python pset 1
# (c) Aaron Bauer 2021

import traceback, os, sys
from inspect import signature

def test_case(test, val, expected, tol=1e-6):
    if type(expected) is tuple:
        if type(val) is not tuple:
            print("FAILED", test, "expected multiple return values but got", val)
            return False
        if len(val) != len(expected):
            print("FAILED", test, "expected", len(expected), "return values but got", len(val))
        for v, e in zip(val, expected):
            if abs(v - e) > tol:
                print("FAILED", test, "expected value was", expected, "but got", val)
                return False
    else:
        if abs(val - expected) > tol:
            print("FAILED", test, "expected value was", expected, "but got", val)
            return False
    return True

if __name__ == "__main__":
    if not os.path.exists("pset1.py"):
        print("ERROR pset1.py not found")
        sys.exit(1)

    # test problem 1: min3, max3, spread3
    try:
        from pset1 import min3, max3, spread3
        try:
            results = []
            results.append(test_case("min3(0, 0, 0)", min3(0, 0, 0), 0))
            results.append(test_case("min3(1, 0, 0)", min3(1, 0, 0), 0))
            results.append(test_case("min3(0, 1, 0)", min3(0, 1, 0), 0))
            results.append(test_case("min3(0, 0, 1)", min3(0, 0, 1), 0))
            results.append(test_case("min3(1, 1, 0)", min3(1, 1, 0), 0))
            results.append(test_case("min3(1, 0, 1)", min3(1, 0, 1), 0))
            results.append(test_case("min3(0, 1, 1)", min3(0, 1, 1), 0))
            results.append(test_case("min3(1, -1, 0)", min3(1, -1, 0), -1))
            if all(results):
                print("PASSED min3")
            else:
                print("FAILED {} of {} min3 tests".format(len([x for x in results if not x]), len(results)))
        except:
            print("ERROR in min3")
            traceback.print_exc()

        try:
            results = []
            results.append(test_case("max3(0, 0, 0)", max3(0, 0, 0), 0))
            results.append(test_case("max3(1, 0, 0)", max3(1, 0, 0), 1))
            results.append(test_case("max3(0, 1, 0)", max3(0, 1, 0), 1))
            results.append(test_case("max3(0, 0, 1)", max3(0, 0, 1), 1))
            results.append(test_case("max3(1, 1, 0)", max3(1, 1, 0), 1))
            results.append(test_case("max3(1, 0, 1)", max3(1, 0, 1), 1))
            results.append(test_case("max3(0, 1, 1)", max3(0, 1, 1), 1))
            results.append(test_case("max3(-1, 1, 0)", max3(-1, 1, 0), 1))
            if all(results):
                print("PASSED max3")
            else:
                print("FAILED {} of {} max3 tests".format(len([x for x in results if not x]), len(results)))
        except:
            print("ERROR in max3")
            traceback.print_exc()
        try:
            results = []
            results.append(test_case("spread3(0, 0, 0)", spread3(0, 0, 0), 0))
            results.append(test_case("spread3(1, -1, 0)", spread3(1, -1, 0), 2))
            results.append(test_case("spread3(-1, 1, 0)", spread3(-1, 1, 0), 2))
            results.append(test_case("spread3(-1, 0, 1)", spread3(-1, 0, 1), 2))
            results.append(test_case("spread3(0, 1, -1)", spread3(0, 1, -1), 2))
            if all(results):
                print("PASSED spread3")
            else:
                print("FAILED {} of {} spread3 tests".format(len([x for x in results if not x]), len(results)))
        except:
            print("ERROR in spread3")
            traceback.print_exc()
    except:
        print("ERROR when importing min3, max3, and spread3")
        traceback.print_exc()

    # test problem 2: convert_angle
    try:
        from pset1 import convert_angle
        if len(signature(convert_angle).parameters) == 3:
            print("found 3-parameter version of convert_angle")
            try:
                results = []
                results.append(test_case("convert_angle(90, 6, 36)", convert_angle(90, 6, 36), 90.11))
                results.append(test_case("convert_angle(-90, 6, 36)", convert_angle(-90, 6, 36), -90.11))
                results.append(test_case("convert_angle(-0.0, 30, 45)", convert_angle(-0.0, 30, 45), -0.5125))
                if all(results):
                    print("convert_angle OK, move on to part b")
                else:
                    print("FAILED {} of {} convert_angle tests".format(len([x for x in results if not x]), len(results)))
            except:
                print("ERROR in convert_angle")
                traceback.print_exc()
        elif len(signature(convert_angle).parameters) == 4:
            print("found 4-parameter version of convert_angle")
            try:
                results = []
                results.append(test_case("convert_angle(90, 6, 36, True)", convert_angle(90, 6, 36, True), 1.57271618897))
                results.append(test_case("convert_angle(-90, 6, 36, True)", convert_angle(-90, 6, 36, True), -1.57271618897))
                results.append(test_case("convert_angle(90, 6, 36, False)", convert_angle(90, 6, 36, False), 90.11))
                results.append(test_case("convert_angle(-90, 6, 36, False)", convert_angle(-90, 6, 36, False), -90.11))
                if all(results):
                    print("convert_angle OK, move on to part c")
                else:
                    print("FAILED {} of {} convert_angle tests".format(len([x for x in results if not x]), len(results)))
            except:
                print("ERROR in convert_angle")
                traceback.print_exc()
        elif len(signature(convert_angle).parameters) == 5:
            print("found 5-parameter version of convert_angle")
            try:
                results = []
                results.append(test_case("convert_angle(90, 6, 36, False, False)", convert_angle(90, 6, 36, False, False), 90.11))
                results.append(test_case("convert_angle(90, 6, 36, True, False)", convert_angle(90, 6, 36, True, False), 1.57271618897))
                results.append(test_case("convert_angle(90, 6, 36, False, True)", convert_angle(90, 6, 36, False, True), 90.11))
                results.append(test_case("convert_angle(90, 6, 36, True, True)", convert_angle(90, 6, 36, True, True), 1.57271618897))
                results.append(test_case("convert_angle(-90, 6, 36, False, False)", convert_angle(-90, 6, 36, False, False), -90.11))
                results.append(test_case("convert_angle(-90, 6, 36, True, False)", convert_angle(-90, 6, 36, True, False), -1.57271618897))
                results.append(test_case("convert_angle(-90, 6, 36, False, True)", convert_angle(-90, 6, 36, False, True), 269.89))
                results.append(test_case("convert_angle(-90, 6, 36, True, True)", convert_angle(-90, 6, 36, True, True), 4.71046911821))
                results.append(test_case("convert_angle(540, 0, 0, False, True)", convert_angle(540, 0, 0, False, True), 180.0))
                results.append(test_case("convert_angle(-0.0, 30, 45, False, False)", convert_angle(-0.0, 30, 45, False, False), -0.5125))
                if all(results):
                    print("PASSED convert_angle")
                else:
                    print("FAILED {} of {} convert_angle tests".format(len([x for x in results if not x]), len(results)))
            except:
                print("ERROR in convert_angle")
                traceback.print_exc()
        else:
            print("ERROR convert_angle takes unexpected number of parameters:", len(signature(convert_angle).parameters))
    except:
        print("ERROR when importing convert_angle")
        traceback.print_exc()

    # test problem 3: RADec_to_AltAz, great_circle_dist
    try:
        from pset1 import RADec_to_AltAz, great_circle_dist
        try:
            results = []
            results.append(test_case("RADec_to_AltAz(156.57025,  24.94131, 253.08608, 34.0727, 2019, 6, 12, 5)", 
                                     RADec_to_AltAz(156.57025,  24.94131, 253.08608, 34.0727, 2019, 6, 12, 5), (28.123690192220057, 282.4588963549189)))
            results.append(test_case("RADec_to_AltAz(238.40339, -19.19939, 253.08608, 34.0727, 2019, 7, 12, 5)", 
                                     RADec_to_AltAz(238.40339, -19.19939, 253.08608, 34.0727, 2019, 7, 12, 5), (33.58536145586729, 202.2271876968573)))
            if all(results):
                print("PASSED RADec_to_AltAz")
            else:
                print("FAILED {} of {} RADec_to_AltAz tests".format(len([x for x in results if not x]), len(results)))
        except:
            print("ERROR in RADec_to_AltAz")
            traceback.print_exc()
        try:
            results = []
            # NMT to CUB; Beijing to DC; Rome to Buenos Aires
            results.append(test_case("great_circle_dist(34.07, -106.9, 40.01, -105.3)", great_circle_dist(34.07, -106.9, 40.01, -105.3), 675555.6199438482, tol=1))
            results.append(test_case("great_circle_dist(39.95, 116.4, 38.92, -77.02)", great_circle_dist(39.95, 116.4, 38.92, -77.02), 11139564.73379008, tol=1))
            results.append(test_case("great_circle_dist(41.91, 12.5, -34.61, -58.45)", great_circle_dist(41.91, 12.5, -34.61, -58.45), 11157229.894962717, tol=1))
            if all(results):
                print("PASSED great_circle_dist")
            else:
                print("FAILED {} of {} great_circle_dist tests".format(len([x for x in results if not x]), len(results)))
        except:
            print("ERROR in great_circle_dist")
            traceback.print_exc()
    except:
        print("ERROR when importing RADec_to_AltAz and great_circle_dist")
        traceback.print_exc()
