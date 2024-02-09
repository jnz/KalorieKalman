KalorieKalman
=============

![LOGO](img/kaloriekalman_small.png)

KalorieKalman is a script package designed to empower individuals on their
fitness journey by providing estimations of total daily energy
expenditure (TDEE), tracking weight variations, and monitoring daily calorie
deficits. Utilizing the robustness of a Kalman filter, KalorieKalman offers an
innovative approach to health and fitness tracking, making it easier for users
to achieve their dietary goals.

Features
--------

    TDEE Estimation: Utilize advanced algorithms to accurately estimate your daily energy expenditure, taking into account various factors for improved accuracy.
    Weight Tracking: Monitor your weight changes over time, providing insights into your progress and helping you stay on track.
    Calorie Deficit Calculation: Automatically calculate your daily calorie deficit, enabling you to understand the effectiveness of your diet and exercise plans.
    Kalman Filter Precision: Incorporate the Kalman filter for enhanced accuracy in predictions and tracking, ensuring you have reliable data to base your decisions on.

Getting Started
===============

Prerequisites
-------------

Before you begin, ensure you have the following installed:

* MATLAB 2017 or higher

Installation (MATLAB)
---------------------

Clone the repository to your local machine:

    git clone https://github.com/yourusername/KalorieKalman.git

Navigate to the cloned repository in MATLAB:

    cd KalorieKalman/matlab

Collect data from iPhone / Apple Health
---------------------------------------

Install the app from the Apple App Store:

    https://apps.apple.com/us/app/simple-health-export-csv/id1535380115

Export the following data types as `.csv` files:

    Body -> BodyMass

    Nutrition -> DietaryEnergyConsumed

Note: this requires tracking your consumed calories (for example with
MyFitnessPal) and body mass (for example with a Bluetooth enabled scale).

Place the `.csv` files in the `matlab` folder and run in MATLAB the `estimateTDEE.m` script.

Contributing
------------

Contributions to KalorieKalman are welcome! If you have suggestions for
improvements or bug fixes, please feel free to fork the repository and submit a
pull request.
