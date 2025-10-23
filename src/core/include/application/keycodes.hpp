#pragma once
#include "exceptions.hpp"

enum Key {
    KEY_A = 65, KEY_B = 66, KEY_C = 67, KEY_D = 68, KEY_E = 69, KEY_F = 70,
    KEY_G = 71, KEY_H = 72, KEY_I = 73, KEY_J = 74, KEY_K = 75, KEY_L = 76,
    KEY_M = 77, KEY_N = 78, KEY_O = 79, KEY_P = 80, KEY_Q = 81, KEY_R = 82,
    KEY_S = 83, KEY_T = 84, KEY_U = 85, KEY_V = 86, KEY_W = 87, KEY_X = 88,
    KEY_Y = 89, KEY_Z = 90,

    KEY_0 = 48, KEY_1 = 49, KEY_2 = 50, KEY_3 = 51, KEY_4 = 52,
    KEY_5 = 53, KEY_6 = 54, KEY_7 = 55, KEY_8 = 56, KEY_9 = 57,

    KEY_SPACE = 32,
    KEY_ESCAPE = 256,
    KEY_ENTER = 257,
    KEY_TAB = 258,
    KEY_BACKSPACE = 259,

    KEY_ARROW_RIGHT = 262,
    KEY_ARROW_LEFT = 263,
    KEY_ARROW_DOWN = 264,
    KEY_ARROW_UP = 265
};

enum MouseButton {
	MOUSE_BUTTON_LEFT = 0,
	MOUSE_BUTTON_RIGHT = 1,
	MOUSE_BUTTON_MIDDLE = 2
};

inline string mouseButtonToString(MouseButton button) {
	switch (button) {
		case MOUSE_BUTTON_LEFT: return "LEFT_BUTTON";
		case MOUSE_BUTTON_RIGHT: return "RIGHT_BUTTON";
		case MOUSE_BUTTON_MIDDLE: return "MIDDLE_BUTTON";
		default: return "UNKNOWN_BUTTON";
	}
}


inline string keyToString(Key key) {
    switch (key) {
        case KEY_A: return "A";
        case KEY_B: return "B";
        case KEY_C: return "C";
        case KEY_D: return "D";
        case KEY_E: return "E";
        case KEY_F: return "F";
        case KEY_G: return "G";
        case KEY_H: return "H";
        case KEY_I: return "I";
        case KEY_J: return "J";
        case KEY_K: return "K";
        case KEY_L: return "L";
        case KEY_M: return "M";
        case KEY_N: return "N";
        case KEY_O: return "O";
        case KEY_P: return "P";
        case KEY_Q: return "Q";
        case KEY_R: return "R";
        case KEY_S: return "S";
        case KEY_T: return "T";
        case KEY_U: return "U";
        case KEY_V: return "V";
        case KEY_W: return "W";
        case KEY_X: return "X";
        case KEY_Y: return "Y";
        case KEY_Z: return "Z";
        case KEY_SPACE: return "SPACE";
        case KEY_ESCAPE: return "ESCAPE";
        case KEY_ENTER: return "ENTER";
        case KEY_TAB: return "TAB";
        case KEY_BACKSPACE: return "BACKSPACE";
        case KEY_ARROW_RIGHT: return "RIGHT";
        case KEY_ARROW_LEFT: return "LEFT";
        case KEY_ARROW_DOWN: return "DOWN";
        case KEY_ARROW_UP: return "UP";
        default: return "UNKNOWN_KEY";
    }
}
