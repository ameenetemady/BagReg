# BagReg
## Introduction
+ BagReg is a novel learning-based protein inference algorithm.
+ BagReg can fully utilize the latent information in the input data.
+ BagReg is competitive with the state-of-the-art protein inference algorithms.

BagReg can be divided into three major phases: feature extraction, prediction model construction and prediction result combination. The method firstly artificially extracts five features from the input data, and then chooses each feature as the class feature to separately build models to predict the presence probabilities of proteins. Finally, the weak results from five prediction models are aggregated to obtain the final result.
## Installation
Our algorithm is writen in Java.
## Usage
Import the code into NetBeans IDE, and then add some datasets, by the way, you should change the path in RAAlgorithm.java. After that, complie and run the code. Finally, the result will be obtained.

Note that the classification method could be manually changed.
