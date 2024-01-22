//todo: continue here
#include <array>
#include "custom_datatypes.h"
//this function calculates the energy need for driving along a certain path
//void getEnergy_v2(SMAsol_struct* sol, double b_model[6], int sol_size[2], model_struct* model, soldata_struct soldata)
double getEnergy_v2(SMAsol_struct* sol, double totalpathlength) //Dummy
{
    return totalpathlength * 5; //Dummy
//  coder::array<double, 2U> b_L_Acc_data;
//  coder::array<double, 2U> b_g_lengthC_data;
//  coder::array<double, 2U> b_g_lengthL_data;
//  coder::array<double, 2U> b_model_h_b;
//  coder::array<double, 2U> c_L_Acc_data;
//  double E_Acc_data[5];
//  double E_per_deg_data[5];
//  double L_Acc_data[5];
//  double a_Acc_data[5];
//  double v_Max_during_Acc_data[5];
//  double t[2];
//  double E_total;
//  double T;
//  double a_tmp;
//  double alphaAll_idx_0;
//  double alphaAll_idx_1;
//  double alphaAll_idx_2;
//  double alphaAll_idx_3;
//  double c_circ_idx_0;
//  double c_circ_idx_1;
//  double c_straight_idx_0;
//  double c_straight_idx_1;
//  double d;
//  double discriminant_tmp;
//  double segmentcounter;
//  double v;
//  double v_allowed;
//  double v_start;
//  int E_per_deg_size[2];
//  int L_Acc_tmp;
//  int b_loop_ub;
//  int i;
//  int i1;
//  int i2;
//  int loop_ub;
//  //  variable initialization
//  E_total = 0.0;
//  //  total energy consumption
//  v = 0.0;
//  //  velocity at given time t
//  v_allowed = model_h_s_Con_v_Max * 0.95;
//  //  maximum allowed velocity
//  segmentcounter = 1.0;
//  //  counter straight lines
//  c_straight_idx_0 = 1.0;
//  c_circ_idx_0 = 1.0;
//  //  counter circles
//  // Initialize for Matlab Coder
//  //     %% Energy Demand during Acceleration
//  // model parameters are interpolated for correct value
//  coder::interp1(model_h_R.d, model_h_b.ED, g_r_data, g_r_size, E_per_deg_data,
//                 E_per_deg_size);
//  //  energy consumption per degree
//  coder::interp1(model_h_R.d, model_h_Acc.E, g_r_data, g_r_size, E_Acc_data,
//                 E_per_deg_size);
//  //  energy cons. acceleration
//  coder::interp1(model_h_R.d, model_h_Acc.L, g_r_data, g_r_size, L_Acc_data,
//                 E_per_deg_size);
//  //  distance ''
//  coder::interp1(model_h_R.d, model_h_Acc.a, g_r_data, g_r_size, a_Acc_data,
//                 E_per_deg_size);
//  //  amount of ''
//  coder::interp1(model_h_R.d, model_h_Acc.v_Max, g_r_data, g_r_size,
//                 v_Max_during_Acc_data, E_per_deg_size);
//  //  max. velocity during ''
//  //  acc. in a straight line
//  //  vmax ''
//  //  distance ''
//  //  energy consumption ''
//  T = 0.0;
//  //  Ampere over time during acceleration
//  //  loop until acceleration done
//  while (v < v_allowed) {
//    //  differentiate straight line, circular line
//    if (rt_remd_snf(segmentcounter, 2.0) != 0.0) {
//      // straight
//      // segment length
//      d = g_lengthL_data[static_cast<int>(c_straight_idx_0) - 1];
//      if (std::abs(d) < 1.0E-5) {
//        //  if segment too short, it is skipped over
//        c_straight_idx_0++;
//      } else {
//        v_start = v;
//        // initial velocity
//        //  Calculate the coefficients of the quadratic equation
//        alphaAll_idx_0 = 0.5 * model_h_s_Acc_a;
//        //  Calculate the discriminant
//        alphaAll_idx_2 = v * v - 4.0 * alphaAll_idx_0 * -d;
//        //  Calculate the roots
//        if (alphaAll_idx_2 >= 0.0) {
//          alphaAll_idx_1 =
//              (-v + std::sqrt(alphaAll_idx_2)) / (2.0 * alphaAll_idx_0);
//        } else {
//          //  No real roots, set t to 0
//          alphaAll_idx_1 = 0.0;
//        }
//        //  Calculate the final velocity
//        v += model_h_s_Acc_a * alphaAll_idx_1;
//        if (v < model_h_s_Con_v_Max) {
//          //  check if value of v allowed
//          // This function calculates E, T and increments E_total
//          // Do not use in cases where formulae are different!
//          E_total += d / model_h_s_Acc_L * model_h_s_Acc_E;
//          T += alphaAll_idx_1;
//          c_straight_idx_0++;
//          // increment counter
//        } else {
//          //  if v too high
//          alphaAll_idx_1 = (model_h_s_Con_v_Max - v_start) / model_h_s_Acc_a;
//          // correct t,s,v
//          alphaAll_idx_2 =
//              v_start * alphaAll_idx_1 +
//              0.5 * model_h_s_Acc_a * (alphaAll_idx_1 * alphaAll_idx_1);
//          v = model_h_s_Con_v_Max;
//          // This function calculates E, T and increments E_total
//          // Do not use in cases where formulae are different!
//          E_total += alphaAll_idx_2 / model_h_s_Acc_L * model_h_s_Acc_E;
//          T += alphaAll_idx_1;
//          //  The distance value is adapted. Oftentimes, there are
//          //  long segments with acceleration and deceleration
//          g_lengthL_data[static_cast<int>(c_straight_idx_0) - 1] =
//              d - alphaAll_idx_2;
//          // adapt segment length
//        }
//      }
//      //  Circles. Steps are analogous to straight line
//    } else {
//      d = g_degC_data[static_cast<int>(c_circ_idx_0) - 1];
//      if (!(std::abs(d) < 1.0E-10)) {
//        v_start = v;
//        //  Calculate the coefficients of the quadratic equation
//        alphaAll_idx_0 = a_Acc_data[static_cast<int>(c_circ_idx_0) - 1];
//        a_tmp = 0.5 * alphaAll_idx_0;
//        //  Calculate the discriminant
//        discriminant_tmp = g_lengthC_data[static_cast<int>(c_circ_idx_0) - 1];
//        alphaAll_idx_2 = v * v - 4.0 * a_tmp * -discriminant_tmp;
//        //  Calculate the roots
//        if (alphaAll_idx_2 >= 0.0) {
//          alphaAll_idx_1 = (-v + std::sqrt(alphaAll_idx_2)) / (2.0 * a_tmp);
//        } else {
//          //  No real roots, set t to 0
//          alphaAll_idx_1 = 0.0;
//        }
//        //  Calculate the final velocity
//        v += alphaAll_idx_0 * alphaAll_idx_1;
//        alphaAll_idx_3 =
//            v_Max_during_Acc_data[static_cast<int>(c_circ_idx_0) - 1];
//        if (v < alphaAll_idx_3) {
//          //  check if v allowed
//          // This function calculates E, T and increments E_total
//          // Do not use in cases where formulae are different!
//          E_total += discriminant_tmp /
//                     L_Acc_data[static_cast<int>(c_circ_idx_0) - 1] *
//                     E_Acc_data[static_cast<int>(c_circ_idx_0) - 1];
//          T += alphaAll_idx_1;
//        } else if (v_start > alphaAll_idx_3) {
//          // if vstart too high
//          v = alphaAll_idx_3;
//          T += discriminant_tmp / alphaAll_idx_3;
//          E_total += E_per_deg_data[static_cast<int>(c_circ_idx_0) - 1] * d;
//        } else {
//          //  if v too high
//          v = alphaAll_idx_3;
//          // correct v,s,t
//          alphaAll_idx_1 = (alphaAll_idx_3 - v_start) / alphaAll_idx_0;
//          // time needed for acc.
//          alphaAll_idx_2 = v_start * alphaAll_idx_1 +
//                           a_tmp * (alphaAll_idx_1 * alphaAll_idx_1);
//          //  time needed to drive with constant velocity
//          // angle corresponding to s
//          // E consumption during acc.
//          // E consumption during deceleration
//          alphaAll_idx_0 = discriminant_tmp - alphaAll_idx_2;
//          E_total =
//              (E_total + alphaAll_idx_2 /
//                             L_Acc_data[static_cast<int>(c_circ_idx_0) - 1] *
//                             E_Acc_data[static_cast<int>(c_circ_idx_0) - 1]) +
//              E_per_deg_data[static_cast<int>(c_circ_idx_0) - 1] *
//                  (alphaAll_idx_0 / discriminant_tmp * d);
//          T = (T + alphaAll_idx_1) + alphaAll_idx_0 / alphaAll_idx_3;
//        }
//      }
//      c_circ_idx_0++;
//    }
//    segmentcounter++;
//  }
//  sol2.L1 = E_total;
//  sol2.T1 = T;
//  //     %% Energy consumption during deceleration
//  v = 0.0;
//  segmentcounter = 1.0;
//  c_straight_idx_1 = g_lengthL_size[1];
//  c_circ_idx_1 = g_lengthC_size[1];
//  coder::interp1(model_h_R.d, model_h_Slw.E, g_r_data, g_r_size, E_Acc_data,
//                 E_per_deg_size);
//  coder::interp1(model_h_R.d, model_h_Slw.L, g_r_data, g_r_size, L_Acc_data,
//                 E_per_deg_size);
//  coder::interp1(model_h_R.d, model_h_Slw.a, g_r_data, g_r_size, a_Acc_data,
//                 E_per_deg_size);
//  //  Ampere over time for deceleration
//  while (v < v_allowed) {
//    //  straight line
//    if (rt_remd_snf(segmentcounter, 2.0) != 0.0) {
//      //  segment too short, skipped over
//      d = g_lengthL_data[static_cast<int>(c_straight_idx_1) - 1];
//      if (std::abs(d) < 1.0E-10) {
//        c_straight_idx_1--;
//      } else {
//        v_start = v;
//        //  Calculate the coefficients of the quadratic equation
//        alphaAll_idx_0 = 0.5 * model_h_s_Slw_a;
//        //  Calculate the discriminant
//        alphaAll_idx_2 = v * v - 4.0 * alphaAll_idx_0 * -d;
//        //  Calculate the roots
//        if (alphaAll_idx_2 >= 0.0) {
//          alphaAll_idx_1 =
//              (-v + std::sqrt(alphaAll_idx_2)) / (2.0 * alphaAll_idx_0);
//        } else {
//          //  No real roots, set t to 0
//          alphaAll_idx_1 = 0.0;
//        }
//        //  Calculate the final velocity
//        v += model_h_s_Slw_a * alphaAll_idx_1;
//        if (v < model_h_s_Con_v_Max) {
//          //  check velocity
//          // This function calculates E, T and increments E_total
//          // Do not use in cases where formulae are different!
//          E_total += d / model_h_s_Slw_L * model_h_s_Slw_E;
//          T += alphaAll_idx_1;
//          c_straight_idx_1--;
//        } else {
//          //  if velocity too high
//          alphaAll_idx_1 = (model_h_s_Con_v_Max - v_start) / model_h_s_Slw_a;
//          alphaAll_idx_2 =
//              v_start * alphaAll_idx_1 +
//              0.5 * model_h_s_Slw_a * (alphaAll_idx_1 * alphaAll_idx_1);
//          v = model_h_s_Con_v_Max;
//          // This function calculates E, T and increments E_total
//          // Do not use in cases where formulae are different!
//          E_total += alphaAll_idx_2 / model_h_s_Slw_L * model_h_s_Slw_E;
//          T += alphaAll_idx_1;
//          g_lengthL_data[static_cast<int>(c_straight_idx_1) - 1] =
//              d - alphaAll_idx_2;
//        }
//      }
//      //  circles
//    } else {
//      d = g_degC_data[static_cast<int>(c_circ_idx_1) - 1];
//      if (!(std::abs(d) < 1.0E-10)) {
//        v_start = v;
//        //  Calculate the coefficients of the quadratic equation
//        alphaAll_idx_0 = a_Acc_data[static_cast<int>(c_circ_idx_1) - 1];
//        a_tmp = 0.5 * alphaAll_idx_0;
//        //  Calculate the discriminant
//        discriminant_tmp = g_lengthC_data[static_cast<int>(c_circ_idx_1) - 1];
//        alphaAll_idx_2 = v * v - 4.0 * a_tmp * -discriminant_tmp;
//        //  Calculate the roots
//        if (alphaAll_idx_2 >= 0.0) {
//          alphaAll_idx_1 = (-v + std::sqrt(alphaAll_idx_2)) / (2.0 * a_tmp);
//        } else {
//          //  No real roots, set t to 0
//          alphaAll_idx_1 = 0.0;
//        }
//        //  Calculate the final velocity
//        v += alphaAll_idx_0 * alphaAll_idx_1;
//        if (v < model_h_s_Con_v_Max) {
//          //  check velocity
//          E_total += discriminant_tmp /
//                     L_Acc_data[static_cast<int>(c_circ_idx_1) - 1] *
//                     E_Acc_data[static_cast<int>(c_circ_idx_1) - 1];
//          T += alphaAll_idx_1;
//        } else {
//          //  if velocity too high
//          v = model_h_s_Con_v_Max;
//          alphaAll_idx_1 = (model_h_s_Con_v_Max - v_start) / alphaAll_idx_0;
//          alphaAll_idx_2 = v_start * alphaAll_idx_1 +
//                           a_tmp * (alphaAll_idx_1 * alphaAll_idx_1);
//          // energy consumption during acc.
//          //  energy cons. while driving with constant velocity
//          alphaAll_idx_0 = discriminant_tmp - alphaAll_idx_2;
//          E_total =
//              (E_total + alphaAll_idx_2 /
//                             L_Acc_data[static_cast<int>(c_circ_idx_1) - 1] *
//                             E_Acc_data[static_cast<int>(c_circ_idx_1) - 1]) +
//              E_per_deg_data[static_cast<int>(c_circ_idx_1) - 1] *
//                  (alphaAll_idx_0 / discriminant_tmp * d);
//          T = (T + alphaAll_idx_1) + alphaAll_idx_0 / model_h_s_Con_v_Max;
//        }
//      }
//      c_circ_idx_1--;
//    }
//    segmentcounter++;
//  }
//  sol2.L2 = E_total - sol2.L1;
//  // sol2.T2 = T - sol2.T1;
//  //  Energy need: Driving with constant velocity
//  //  previous energy need + circular segments + straight segments
//  if (c_circ_idx_0 > c_circ_idx_1) {
//    i = 0;
//    i1 = 0;
//  } else {
//    i = static_cast<int>(c_circ_idx_0) - 1;
//    i1 = static_cast<int>(c_circ_idx_1);
//  }
//  if (c_straight_idx_0 > c_straight_idx_1) {
//    i2 = 0;
//    loop_ub = 0;
//  } else {
//    i2 = static_cast<int>(c_straight_idx_0) - 1;
//    loop_ub = static_cast<int>(c_straight_idx_1);
//  }
//  b_loop_ub = i1 - i;
//  for (i1 = 0; i1 < b_loop_ub; i1++) {
//    L_Acc_tmp = i + i1;
//    L_Acc_data[i1] = E_per_deg_data[L_Acc_tmp] * g_degC_data[L_Acc_tmp];
//  }
//  loop_ub -= i2;
//  b_model_h_b.set_size(1, loop_ub);
//  for (i = 0; i < loop_ub; i++) {
//    b_model_h_b[i] = model_h_b.ELs * g_lengthL_data[i2 + i];
//  }
//  b_L_Acc_data.set(&L_Acc_data[0], 1, b_loop_ub);
//  sol2.L = (E_total + coder::sum(b_L_Acc_data)) + coder::sum(b_model_h_b);
//  //  turn is performed while standing.
//  t[0] = g_startTurn;
//  t[1] = g_endTurn;
//  if (coder::any(t)) {
//    double dv[2];
//    // checks if there is standing rotation
//    alphaAll_idx_0 = sol2.XS[1] - sol2.XS[0];
//    alphaAll_idx_2 = sol2.YS[1] - sol2.YS[0];
//    alphaAll_idx_1 = sol2.XS[8] - sol2.XS[9];
//    alphaAll_idx_3 = sol2.YS[8] - sol2.YS[9];
//    if (g_startTurn != 0.0) {
//      //  if amount of standing rotation is not zero
//      // This function calculates the energy consumption during a standing
//      // rotation try get angle in relation to x-axis copied from
//      // ParseSolution_SMA, because Matlab won't find it otherwise
//      t[0] = alphaAll_idx_0;
//      dv[0] = 1.0;
//      t[1] = alphaAll_idx_2;
//      dv[1] = 0.0;
//      if (std::isnan(alphaAll_idx_2)) {
//        d = rtNaN;
//      } else if (alphaAll_idx_2 < 0.0) {
//        d = -1.0;
//      } else {
//        d = (alphaAll_idx_2 > 0.0);
//      }
//      alphaAll_idx_0 = std::acos((alphaAll_idx_0 + alphaAll_idx_2 * 0.0) /
//                                 (coder::b_norm(t) * coder::b_norm(dv))) *
//                       d * 180.0 / 3.1415926535897931;
//      // radtodeg conversion
//      // get angular difference between beta (desired) and alpha (real)
//      // catch
//      // end
//      // handle cases with alpha<0 or alpha >180
//      if (alphaAll_idx_0 < 0.0) {
//        alphaAll_idx_0 += 360.0;
//      }
//      alphaAll_idx_0 = std::abs(model_SA - alphaAll_idx_0);
//      if (alphaAll_idx_0 > 180.0) {
//        alphaAll_idx_0 = 360.0 - alphaAll_idx_0;
//      }
//      alphaAll_idx_2 = model_h_c0_Acc.D + model_h_c0_Slw.D;
//      // calculate energy need
//      if (alphaAll_idx_0 > alphaAll_idx_2) {
//        alphaAll_idx_0 = (alphaAll_idx_0 - model_h_c0_Acc.D) - model_h_c0_Slw.D;
//        alphaAll_idx_0 = (model_h_c0_Acc.E + model_h_c0_Slw.E) +
//                         alphaAll_idx_0 * model_h_c0_Con.b_ED;
//      } else {
//        alphaAll_idx_0 = alphaAll_idx_0 * (model_h_c0_Acc.D / alphaAll_idx_2) /
//                             model_h_c0_Acc.D * model_h_c0_Acc.E +
//                         alphaAll_idx_0 * (model_h_c0_Slw.D / alphaAll_idx_2) /
//                             model_h_c0_Slw.D * model_h_c0_Slw.E;
//      }
//      // add result to L
//      //  Add result to L
//      sol2.L1 += alphaAll_idx_0;
//      sol2.L += alphaAll_idx_0;
//    }
//    if (g_endTurn != 0.0) {
//      //  if amount of standing rotation is not zero
//      // This function calculates the energy consumption during a standing
//      // rotation try get angle in relation to x-axis copied from
//      // ParseSolution_SMA, because Matlab won't find it otherwise
//      t[0] = alphaAll_idx_1;
//      dv[0] = 1.0;
//      t[1] = alphaAll_idx_3;
//      dv[1] = 0.0;
//      if (std::isnan(alphaAll_idx_3)) {
//        d = rtNaN;
//      } else if (alphaAll_idx_3 < 0.0) {
//        d = -1.0;
//      } else {
//        d = (alphaAll_idx_3 > 0.0);
//      }
//      alphaAll_idx_0 = std::acos((alphaAll_idx_1 + alphaAll_idx_3 * 0.0) /
//                                 (coder::b_norm(t) * coder::b_norm(dv))) *
//                       d * 180.0 / 3.1415926535897931;
//      // radtodeg conversion
//      // get angular difference between beta (desired) and alpha (real)
//      // catch
//      // end
//      // handle cases with alpha<0 or alpha >180
//      if (alphaAll_idx_0 < 0.0) {
//        alphaAll_idx_0 += 360.0;
//      }
//      alphaAll_idx_0 = std::abs(model_EA - alphaAll_idx_0);
//      if (alphaAll_idx_0 > 180.0) {
//        alphaAll_idx_0 = 360.0 - alphaAll_idx_0;
//      }
//      alphaAll_idx_2 = model_h_c0_Acc.D + model_h_c0_Slw.D;
//      // calculate energy need
//      if (alphaAll_idx_0 > alphaAll_idx_2) {
//        alphaAll_idx_0 = (alphaAll_idx_0 - model_h_c0_Acc.D) - model_h_c0_Slw.D;
//        alphaAll_idx_0 = (model_h_c0_Acc.E + model_h_c0_Slw.E) +
//                         alphaAll_idx_0 * model_h_c0_Con.b_ED;
//      } else {
//        alphaAll_idx_0 = alphaAll_idx_0 * (model_h_c0_Acc.D / alphaAll_idx_2) /
//                             model_h_c0_Acc.D * model_h_c0_Acc.E +
//                         alphaAll_idx_0 * (model_h_c0_Slw.D / alphaAll_idx_2) /
//                             model_h_c0_Slw.D * model_h_c0_Slw.E;
//      }
//      // add result to L
//      //  Add result to L
//      sol2.L2 += alphaAll_idx_0;
//      sol2.L += alphaAll_idx_0;
//    }
//  }
//  //  other interesting data
//  b_g_lengthL_data.set(&g_lengthL_data[0], g_lengthL_size[0],
//                       g_lengthL_size[1]);
//  b_g_lengthC_data.set((double *)&g_lengthC_data[0], g_lengthC_size[0],
//                       g_lengthC_size[1]);
//  sol2.length = coder::sum(b_g_lengthL_data) + coder::sum(b_g_lengthC_data);
//  //  distance
//  //  Energy need over time
//  coder::interp1(model_h_R.d, model_h_b.TD, g_r_data, g_r_size, E_Acc_data,
//                 E_per_deg_size);
//  if (c_circ_idx_0 > c_circ_idx_1) {
//    i = 0;
//    i1 = 0;
//  } else {
//    i = static_cast<int>(c_circ_idx_0) - 1;
//    i1 = static_cast<int>(c_circ_idx_1);
//  }
//  if (c_straight_idx_0 > c_straight_idx_1) {
//    i2 = 0;
//    loop_ub = 0;
//  } else {
//    i2 = static_cast<int>(c_straight_idx_0) - 1;
//    loop_ub = static_cast<int>(c_straight_idx_1);
//  }
//  b_loop_ub = i1 - i;
//  for (i1 = 0; i1 < b_loop_ub; i1++) {
//    L_Acc_tmp = i + i1;
//    L_Acc_data[i1] = E_Acc_data[L_Acc_tmp] * g_degC_data[L_Acc_tmp];
//  }
//  loop_ub -= i2;
//  b_model_h_b.set_size(1, loop_ub);
//  for (i = 0; i < loop_ub; i++) {
//    b_model_h_b[i] = model_h_b.TLs * g_lengthL_data[i2 + i];
//  }
//  c_L_Acc_data.set(&L_Acc_data[0], 1, b_loop_ub);
//  sol2.T = (T + coder::sum(c_L_Acc_data)) + coder::sum(b_model_h_b);
//  //  Computation of Ampere over time for driving with constant velocity as well
//  //  as necessary parameters. docking, undocking from charging station
//  if (model_A_is_CS != 0.0) {
//    sol2.L += model_h_updock.E;
//    sol2.T += model_h_updock.T;
//    sol2.length += model_h_updock.L;
//  }
//  if (model_B_is_CS != 0.0) {
//    sol2.L += model_h_dock.E;
//    sol2.T += model_h_dock.T;
//    sol2.length += model_h_dock.L;
//  }
//  // if eval
//  //     sol2.T_array = T_array;
//  //     sol2.L_array_A = L_array_A;
//  // end
}
//