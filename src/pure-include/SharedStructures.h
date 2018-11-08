/*
 * Copyright (c) 2013 - Idiap Research Institute - http://www.idiap.ch/
 * 
 * See file COPYING for the licence associated to this software
 * 
 * This file is a part of the TopicModels4Video software version: academic-1.0
*/

namespace Shared {

    struct LabeledDetection {
        int x;
        int y;
        int width;
        int height;
        int label;
    };


    // simple rectangles
    struct DetectionBox {
        int x;
        int y;
        int width;
        int height;
    };

    enum LLState {YOUNG, WARNED, ALARMED, REMOVED};
    struct LLInfo {
        LLState previousState;
        LLState currentState;
        int currentFrame;
        DetectionBox bbox;
        int dropFrame;
        int warnFrame;
        int alarmFrame;
    };
}


