## Unscented Kalman Filter Project. [Rubric](https://review.udacity.com/#!/rubrics/783/view) Points
Self-Driving Car Engineer Nanodegree Program

---
### 1. Accuracy for data-1

> For the new version of the project, there is now only one data set "obj_pose-laser-radar-synthetic-input.txt". px, py, vx, vy output coordinates must have an RMSE <= [.09, .10, .40, .30] when using the file: "obj_pose-laser-radar-synthetic-input.txt"

I got

| My Accuracy | Target Accuracy|
|:-------------:|:-------------:|
| 0.0670838     | 0.09 |
| 0.0816007     | 0.10 |
| 0.398697      | 0.40 |
| 0.207379      | 0.30 |

---
### 2. Algorithm-1

> Your Sensor Fusion algorithm follows the general processing flow as taught in the preceding lessons.

I followed the flow taught in the lesson. Most building blocks, including Predict(), AugmentedSigmaPoints(), PredictMeanAndCovariance() , PredictRadarMeasurement() and Update() are actually the codes in the lecture.


---
### 3. Algorithm-2

> Your Kalman Filter algorithm handles the first measurements appropriately.

I handled the first measurements as follows:

1. RADAR:
``` python
float rho = meas_package.raw_measurements_(0);
float phi = meas_package.raw_measurements_(1);
x_ << sin(phi)*rho, cos(phi)*rho, abs(meas_package.raw_measurements_(2)), -phi+M_PI/2, 0;
```

2. LASER
``` python
x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 1, M_PI, 0;
```

Before calling those assiangment, a checking on nonsence sensor data and the type of used sensor are done by:
``` python
// type of used sensor
if (use_radar_==false && (meas_package.sensor_type_ == MeasurementPackage::RADAR))
    return;
if (use_laser_==false && (meas_package.sensor_type_ == MeasurementPackage::LASER))
    return;
// Avoid initializing nonsence data
if (x_[0] < MIN_SENSOR_VALUE && x_[1] < MIN_SENSOR_VALUE)
{
    std::cout << "Avoid initializing nonsence data\n";
    return;
}
```

---
### 4. Algorithm-3

> Your Kalman Filter algorithm first predicts then updates.

Sure it is! The codes structure are below:
``` c++
/*****************************************************************************
 *  Prediction
 ****************************************************************************/
double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
time_us_ = meas_package.timestamp_;

Prediction(dt);

/*****************************************************************************
 *  Update
 ****************************************************************************/

if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);
} else {
    // Laser updates
    UpdateLidar(meas_package);
}
```

---
### 5. Algorithm-4

> Your Kalman Filter can handle radar and lidar measurements.

Sure it is! As shown in the above code snippet, I update with `UpdateRadar()` and `UpdateLidar()` according to what type of sensor is.

---
### 6. Code Efficiency

> Your algorithm should avoid unnecessary calculations.

1. I avoid running the same calculation in `UpdateLidar()`, `UpdateRadar()` and `Prediction()`. More specifically for `Prediction()`, when delta_t is too small (~0), the predicted sigma points can be calucuated more efficiently!
```c++
if (delta_t<0.0001)
{
  // Can just assign values since delta_t is almost 0. So we can save some calculations.
    px_p = p_x;
    py_p = p_y;
    v_p = v;
    yaw_p = yaw;
    yawd_p = yawd;
}
{
  // CTRV state recursion calculation ...
}
//write predicted sigma point into right column
Xsig_pred_(0,i) = px_p;
Xsig_pred_(1,i) = py_p;
Xsig_pred_(2,i) = v_p;
Xsig_pred_(3,i) = yaw_p;
Xsig_pred_(4,i) = yawd_p;
```
2. No complex data structures are used. In fact, the template codes provided are well structured and most variables we need are defined already.
3. No Unnecessary control flow checks.
4. Efficiently compute angle normalization
Using direct calculation instead of the loop method
```c++
if (fabs(x) > M_PI)
{
    x -= round(x / (2.0 * M_PI)) * (2.0 * M_PI);
}
```

