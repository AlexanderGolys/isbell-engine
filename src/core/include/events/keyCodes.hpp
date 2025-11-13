#pragma once

enum class KeyCode {
	UNKNOWN = 0,
	SPACE = 32,
	APOSTROPHE = 39, /* ' */
	COMMA = 44,      /* , */
	MINUS = 45,      /* - */
	PERIOD = 46,     /* . */
	SLASH = 47,      /* / */
	KEY_0 = 48,
	KEY_1 = 49,
	KEY_2 = 50,
	KEY_3 = 51,
	KEY_4 = 52,
	KEY_5 = 53,
	KEY_6 = 54,
	KEY_7 = 55,
	KEY_8 = 56,
	KEY_9 = 57,
	SEMICOLON = 59, /* ; */
	EQUAL = 61,     /* = */
	KEY_A = 65,
	KEY_B = 66,
	KEY_C = 67,
	KEY_D = 68,
	KEY_E = 69,
	KEY_F = 70,
	KEY_G = 71,
	KEY_H = 72,
	KEY_I = 73,
	KEY_J = 74,
	KEY_K = 75,
	KEY_L = 76,
	KEY_M = 77,
	KEY_N = 78,
	KEY_O = 79,
	KEY_P = 80,
	KEY_Q = 81,
	KEY_R = 82,
	KEY_S = 83,
	KEY_T = 84,
	KEY_U = 85,
	KEY_V = 86,
	KEY_W = 87,
	KEY_X = 88,
	KEY_Y = 89,
	KEY_Z = 90,
	LEFT_BRACKET = 91,  /* [ */
	BACKSLASH = 92,     /* \ */
	RIGHT_BRACKET = 93, /* ] */
	GRAVE_ACCENT = 96,  /* ` */
	ESCAPE = 256,
	ENTER = 257,
	TAB = 258,
	BACKSPACE = 259,
	INSERT = 260,
	DELETE_KEY = 261,
	RIGHT_ARROW = 262,
	LEFT_ARROW = 263,
	DOWN_ARROW = 264,
	UP_ARROW = 265,
	PAGE_UP = 266,
	PAGE_DOWN = 267,
	HOME = 268,
	END = 269,
	CAPS_LOCK = 280,
	SCROLL_LOCK = 281,
	NUM_LOCK = 282,
	PRINT_SCREEN = 283
};

enum class MouseButton {
	LEFT = 0,
	RIGHT = 1,
	MIDDLE = 2
};