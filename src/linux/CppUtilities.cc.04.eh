
;; Function void green2redRGB(unsigned char&, unsigned char&, unsigned char&, double, double, double, double) (_Z12green2redRGBRhS_S_dddd)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->3 redirected to 4
Forwarding edge 0->2 to 5 failed.
Forwarding edge 0->2 to 5 failed.
Deleting block 3.
Redirecting jump 31 from 23 to 25.
Edge 5->7 redirected to 8
Forwarding edge 5->6 to 9 failed.
Forwarding edge 5->6 to 9 failed.
Deleting block 7.
Redirecting jump 49 from 23 to 25.
Edge 9->11 redirected to 12
Forwarding edge 9->10 to 16 failed.
Forwarding edge 9->10 to 16 failed.
Deleting block 11.
Edge 12->14 redirected to 15
Forwarding edge 12->13 to 16 failed.
Forwarding edge 12->13 to 16 failed.
Deleting block 14.
Redirecting jump 85 from 23 to 25.
Edge 16->18 redirected to 19
Forwarding edge 16->17 to 25 failed.
Forwarding edge 16->17 to 25 failed.
Redirecting jump 92 from 23 to 25.
Deleting block 18.
Edge 19->21 redirected to 22
Forwarding edge 19->20 to 25 failed.
Forwarding edge 19->20 to 25 failed.
Redirecting jump 99 from 23 to 25.
Deleting block 21.
Merged 22 and 23 without moving.
Merged 22 and 24 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 5 failed.
Forwarding edge 5->6 to 9 failed.
Forwarding edge 9->10 to 16 failed.
Forwarding edge 12->13 to 16 failed.
Forwarding edge 16->17 to 25 failed.
Forwarding edge 19->20 to 25 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 3 failed.
Forwarding edge 3->4 to 6 failed.
Forwarding edge 6->7 to 11 failed.
Forwarding edge 8->9 to 11 failed.
Forwarding edge 11->12 to 16 failed.
Forwarding edge 13->14 to 16 failed.
(note 1 0 15 ("./CppUtilities.cc") 200)

;; Start of basic block 0, registers live: (nil)
(note 15 1 6 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 6 15 7 0 ./CppUtilities.cc:200 (set (reg/v/f:DI 60 [ r ])
        (reg:DI 5 di [ r ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppUtilities.cc:200 (set (reg/v/f:DI 61 [ g ])
        (reg:DI 4 si [ g ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppUtilities.cc:200 (set (reg/v/f:DI 62 [ b ])
        (reg:DI 1 dx [ b ])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppUtilities.cc:200 (set (reg/v:DF 63 [ signal ])
        (reg:DF 21 xmm0 [ signal ])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./CppUtilities.cc:200 (set (reg/v:DF 64 [ min ])
        (reg:DF 22 xmm1 [ min ])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppUtilities.cc:200 (set (reg/v:DF 65 [ max ])
        (reg:DF 23 xmm2 [ max ])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppUtilities.cc:200 (set (reg/v:DF 66 [ med ])
        (reg:DF 24 xmm3 [ med ])) -1 (nil)
    (nil))

(note 13 12 17 0 NOTE_INSN_FUNCTION_BEG)

(note 17 13 18 0 ("./CppUtilities.cc") 204)

(insn 18 17 19 0 ./CppUtilities.cc:204 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 64 [ min ])
            (reg/v:DF 63 [ signal ]))) -1 (nil)
    (nil))

(jump_insn 19 18 131 0 ./CppUtilities.cc:204 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 23)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 131 19 20 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 20 131 21 1 ./CppUtilities.cc:204 (set (pc)
        (label_ref 33)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 21 20 23)

;; Start of basic block 2, registers live: (nil)
(code_label 23 21 24 2 4 "" [1 uses])

(note 24 23 25 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 25 24 26 2 ("./CppUtilities.cc") 206)

(insn 26 25 27 2 ./CppUtilities.cc:206 (set (mem:QI (reg/v/f:DI 60 [ r ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 27 26 28 2 ("./CppUtilities.cc") 207)

(insn 28 27 29 2 ./CppUtilities.cc:207 (set (mem:QI (reg/v/f:DI 61 [ g ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 29 28 30 2 ("./CppUtilities.cc") 208)

(insn 30 29 31 2 ./CppUtilities.cc:208 (set (mem:QI (reg/v/f:DI 62 [ b ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 31 30 32 2 ./CppUtilities.cc:208 (set (pc)
        (label_ref:DI 129)) 548 {jump} (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 32 31 33)

;; Start of basic block 3, registers live: (nil)
(code_label 33 32 34 3 2 "" [1 uses])

(note 34 33 35 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 35 34 36 3 ("./CppUtilities.cc") 210)

(insn 36 35 37 3 ./CppUtilities.cc:210 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 63 [ signal ])
            (reg/v:DF 65 [ max ]))) -1 (nil)
    (nil))

(jump_insn 37 36 133 3 ./CppUtilities.cc:210 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 41)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 133 37 38 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(jump_insn 38 133 39 4 ./CppUtilities.cc:210 (set (pc)
        (label_ref 51)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

(barrier 39 38 41)

;; Start of basic block 5, registers live: (nil)
(code_label 41 39 42 5 8 "" [1 uses])

(note 42 41 43 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 43 42 44 5 ("./CppUtilities.cc") 212)

(insn 44 43 45 5 ./CppUtilities.cc:212 (set (mem:QI (reg/v/f:DI 60 [ r ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 45 44 46 5 ("./CppUtilities.cc") 213)

(insn 46 45 47 5 ./CppUtilities.cc:213 (set (mem:QI (reg/v/f:DI 61 [ g ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 47 46 48 5 ("./CppUtilities.cc") 214)

(insn 48 47 49 5 ./CppUtilities.cc:214 (set (mem:QI (reg/v/f:DI 62 [ b ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 49 48 50 5 ./CppUtilities.cc:214 (set (pc)
        (label_ref:DI 129)) 548 {jump} (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 50 49 51)

;; Start of basic block 6, registers live: (nil)
(code_label 51 50 52 6 6 "" [1 uses])

(note 52 51 53 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 53 52 54 6 ("./CppUtilities.cc") 216)

(insn 54 53 55 6 ./CppUtilities.cc:216 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 63 [ signal ])
            (reg/v:DF 64 [ min ]))) -1 (nil)
    (nil))

(jump_insn 55 54 135 6 ./CppUtilities.cc:216 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 59)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 135 55 56 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(jump_insn 56 135 57 7 ./CppUtilities.cc:216 (set (pc)
        (label_ref 87)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

(barrier 57 56 59)

;; Start of basic block 8, registers live: (nil)
(code_label 59 57 60 8 11 "" [1 uses])

(note 60 59 61 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 61 60 62 8 ./CppUtilities.cc:216 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 66 [ med ])
            (reg/v:DF 63 [ signal ]))) -1 (nil)
    (nil))

(jump_insn 62 61 137 8 ./CppUtilities.cc:216 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 66)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 137 62 63 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(jump_insn 63 137 64 9 ./CppUtilities.cc:216 (set (pc)
        (label_ref 87)) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

(barrier 64 63 66)

;; Start of basic block 10, registers live: (nil)
(code_label 66 64 67 10 13 "" [1 uses])

(note 67 66 68 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 68 67 69 10 ("./CppUtilities.cc") 218)

(insn 69 68 70 10 ./CppUtilities.cc:218 (set (reg:DF 67)
        (minus:DF (reg/v:DF 66 [ med ])
            (reg/v:DF 64 [ min ]))) -1 (nil)
    (nil))

(insn 70 69 71 10 ./CppUtilities.cc:218 (set (reg:DF 69)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 71 70 72 10 ./CppUtilities.cc:218 (set (reg:DF 68)
        (div:DF (reg:DF 69)
            (reg:DF 67))) -1 (nil)
    (nil))

(insn 72 71 73 10 ./CppUtilities.cc:218 (set (reg:DF 70)
        (minus:DF (reg/v:DF 63 [ signal ])
            (reg/v:DF 64 [ min ]))) -1 (nil)
    (nil))

(insn 73 72 74 10 ./CppUtilities.cc:218 (set (reg/v:DF 58 [ t.236 ])
        (mult:DF (reg:DF 68)
            (reg:DF 70))) -1 (nil)
    (nil))

(note 74 73 75 10 ("./CppUtilities.cc") 220)

(insn 75 74 76 10 ./CppUtilities.cc:220 (set (mem:QI (reg/v/f:DI 60 [ r ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 76 75 77 10 ("./CppUtilities.cc") 221)

(insn 77 76 78 10 ./CppUtilities.cc:221 (set (reg:DF 71)
        (mult:DF (reg/v:DF 58 [ t.236 ])
            (reg/v:DF 58 [ t.236 ]))) -1 (nil)
    (nil))

(insn 78 77 79 10 ./CppUtilities.cc:221 (set (reg:DF 72)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 2.55e+2 [0x0.ffp+8])
        (nil)))

(insn 79 78 80 10 ./CppUtilities.cc:221 (set (reg:DF 73)
        (mult:DF (reg:DF 71)
            (reg:DF 72))) -1 (nil)
    (nil))

(insn 80 79 81 10 ./CppUtilities.cc:221 (set (reg:SI 74)
        (fix:SI (reg:DF 73))) -1 (nil)
    (nil))

(insn 81 80 82 10 ./CppUtilities.cc:221 (set (mem:QI (reg/v/f:DI 61 [ g ]) [0 S1 A8])
        (subreg:QI (reg:SI 74) 0)) -1 (nil)
    (nil))

(note 82 81 83 10 ("./CppUtilities.cc") 222)

(insn 83 82 84 10 ./CppUtilities.cc:222 (set (mem:QI (reg/v/f:DI 62 [ b ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 84 83 85 10 ("./CppUtilities.cc") 216)

(jump_insn 85 84 86 10 ./CppUtilities.cc:216 (set (pc)
        (label_ref:DI 129)) 548 {jump} (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 86 85 87)

;; Start of basic block 11, registers live: (nil)
(code_label 87 86 88 11 9 "" [2 uses])

(note 88 87 89 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(note 89 88 90 11 ("./CppUtilities.cc") 224)

(insn 90 89 91 11 ./CppUtilities.cc:224 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 63 [ signal ])
            (reg/v:DF 66 [ med ]))) -1 (nil)
    (nil))

(jump_insn 91 90 139 11 ./CppUtilities.cc:224 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 95)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 6700 [0x1a2c])
        (nil)))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(note 139 91 92 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(jump_insn 92 139 93 12 ./CppUtilities.cc:224 (set (pc)
        (label_ref:DI 129)) 548 {jump} (nil)
    (nil))
;; End of basic block 12, registers live:
 (nil)

(barrier 93 92 95)

;; Start of basic block 13, registers live: (nil)
(code_label 95 93 96 13 15 "" [1 uses])

(note 96 95 97 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 97 96 98 13 ./CppUtilities.cc:224 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 65 [ max ])
            (reg/v:DF 63 [ signal ]))) -1 (nil)
    (nil))

(jump_insn 98 97 141 13 ./CppUtilities.cc:224 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 102)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 6700 [0x1a2c])
        (nil)))
;; End of basic block 13, registers live:
 (nil)

;; Start of basic block 14, registers live: (nil)
(note 141 98 99 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(jump_insn 99 141 100 14 ./CppUtilities.cc:224 (set (pc)
        (label_ref:DI 129)) 548 {jump} (nil)
    (nil))
;; End of basic block 14, registers live:
 (nil)

(barrier 100 99 102)

;; Start of basic block 15, registers live: (nil)
(code_label 102 100 103 15 17 "" [1 uses])

(note 103 102 104 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 104 103 105 15 ("./CppUtilities.cc") 226)

(insn 105 104 106 15 ./CppUtilities.cc:226 (set (reg:DF 75)
        (minus:DF (reg/v:DF 65 [ max ])
            (reg/v:DF 66 [ med ]))) -1 (nil)
    (nil))

(insn 106 105 107 15 ./CppUtilities.cc:226 (set (reg:DF 76)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC2") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF -1.0e+0 [-0x0.8p+1])
        (nil)))

(insn 107 106 108 15 ./CppUtilities.cc:226 (set (reg:DF 77)
        (div:DF (reg:DF 76)
            (reg:DF 75))) -1 (nil)
    (nil))

(insn 108 107 109 15 ./CppUtilities.cc:226 (set (reg:DF 78)
        (minus:DF (reg/v:DF 63 [ signal ])
            (reg/v:DF 66 [ med ]))) -1 (nil)
    (nil))

(insn 109 108 110 15 ./CppUtilities.cc:226 (set (reg:DF 79)
        (mult:DF (reg:DF 77)
            (reg:DF 78))) -1 (nil)
    (nil))

(insn 110 109 111 15 ./CppUtilities.cc:226 (set (reg:DF 80)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 111 110 112 15 ./CppUtilities.cc:226 (set (reg/v:DF 59 [ t ])
        (plus:DF (reg:DF 79)
            (reg:DF 80))) -1 (nil)
    (nil))

(note 112 111 113 15 ("./CppUtilities.cc") 228)

(insn 113 112 114 15 ./CppUtilities.cc:228 (set (reg:DF 81)
        (mult:DF (reg/v:DF 59 [ t ])
            (reg/v:DF 59 [ t ]))) -1 (nil)
    (nil))

(insn 114 113 115 15 ./CppUtilities.cc:228 (set (reg:DF 82)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 2.55e+2 [0x0.ffp+8])
        (nil)))

(insn 115 114 116 15 ./CppUtilities.cc:228 (set (reg:DF 83)
        (mult:DF (reg:DF 81)
            (reg:DF 82))) -1 (nil)
    (nil))

(insn 116 115 117 15 ./CppUtilities.cc:228 (set (reg:SI 84)
        (fix:SI (reg:DF 83))) -1 (nil)
    (nil))

(insn 117 116 118 15 ./CppUtilities.cc:228 (set (mem:QI (reg/v/f:DI 60 [ r ]) [0 S1 A8])
        (subreg:QI (reg:SI 84) 0)) -1 (nil)
    (nil))

(note 118 117 119 15 ("./CppUtilities.cc") 229)

(insn 119 118 120 15 ./CppUtilities.cc:229 (set (mem:QI (reg/v/f:DI 61 [ g ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 120 119 121 15 ("./CppUtilities.cc") 230)

(insn 121 120 124 15 ./CppUtilities.cc:230 (set (mem:QI (reg/v/f:DI 62 [ b ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 15, registers live:
 (nil)

(note 124 121 125 NOTE_INSN_FUNCTION_END)

(note 125 124 129 ("./CppUtilities.cc") 234)

;; Start of basic block 16, registers live: (nil)
(code_label 129 125 144 16 18 "" [5 uses])

(note 144 129 0 16 [bb 16] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 16, registers live:
 (nil)


;; Function void green2redHSV(unsigned char&, unsigned char&, unsigned char&, double, double, double, double) (_Z12green2redHSVRhS_S_dddd)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->3 redirected to 4
Forwarding edge 0->2 to 5 failed.
Forwarding edge 0->2 to 5 failed.
Deleting block 3.
Redirecting jump 28 from 23 to 25.
Edge 5->7 redirected to 8
Forwarding edge 5->6 to 9 failed.
Forwarding edge 5->6 to 9 failed.
Deleting block 7.
Redirecting jump 46 from 23 to 25.
Edge 9->11 redirected to 12
Forwarding edge 9->10 to 16 failed.
Forwarding edge 9->10 to 16 failed.
Deleting block 11.
Edge 12->14 redirected to 15
Forwarding edge 12->13 to 16 failed.
Forwarding edge 12->13 to 16 failed.
Deleting block 14.
Redirecting jump 83 from 23 to 25.
Edge 16->18 redirected to 19
Forwarding edge 16->17 to 25 failed.
Forwarding edge 16->17 to 25 failed.
Redirecting jump 90 from 23 to 25.
Deleting block 18.
Edge 19->21 redirected to 22
Forwarding edge 19->20 to 25 failed.
Forwarding edge 19->20 to 25 failed.
Redirecting jump 97 from 23 to 25.
Deleting block 21.
Merged 22 and 23 without moving.
Merged 22 and 24 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 5 failed.
Forwarding edge 5->6 to 9 failed.
Forwarding edge 9->10 to 16 failed.
Forwarding edge 12->13 to 16 failed.
Forwarding edge 16->17 to 25 failed.
Forwarding edge 19->20 to 25 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 3 failed.
Forwarding edge 3->4 to 6 failed.
Forwarding edge 6->7 to 11 failed.
Forwarding edge 8->9 to 11 failed.
Forwarding edge 11->12 to 16 failed.
Forwarding edge 13->14 to 16 failed.
(note 1 0 12 ("./CppUtilities.cc") 240)

;; Start of basic block 0, registers live: (nil)
(note 12 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 12 4 0 ./CppUtilities.cc:240 (set (reg/v/f:DI 59 [ h ])
        (reg:DI 5 di [ h ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:240 (set (reg/v/f:DI 60 [ s ])
        (reg:DI 4 si [ s ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppUtilities.cc:240 (set (reg/v/f:DI 61 [ v ])
        (reg:DI 1 dx [ v ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppUtilities.cc:240 (set (reg/v:DF 62 [ signal ])
        (reg:DF 21 xmm0 [ signal ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppUtilities.cc:240 (set (reg/v:DF 63 [ min ])
        (reg:DF 22 xmm1 [ min ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppUtilities.cc:240 (set (reg/v:DF 64 [ max ])
        (reg:DF 23 xmm2 [ max ])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppUtilities.cc:240 (set (reg/v:DF 65 [ med ])
        (reg:DF 24 xmm3 [ med ])) -1 (nil)
    (nil))

(note 10 9 14 0 NOTE_INSN_FUNCTION_BEG)

(note 14 10 15 0 ("./CppUtilities.cc") 244)

(insn 15 14 16 0 ./CppUtilities.cc:244 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 63 [ min ])
            (reg/v:DF 62 [ signal ]))) -1 (nil)
    (nil))

(jump_insn 16 15 129 0 ./CppUtilities.cc:244 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 20)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 129 16 17 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 17 129 18 1 ./CppUtilities.cc:244 (set (pc)
        (label_ref 30)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 18 17 20)

;; Start of basic block 2, registers live: (nil)
(code_label 20 18 21 2 25 "" [1 uses])

(note 21 20 22 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 22 21 23 2 ("./CppUtilities.cc") 246)

(insn 23 22 24 2 ./CppUtilities.cc:246 (set (mem:QI (reg/v/f:DI 59 [ h ]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 24 23 25 2 ("./CppUtilities.cc") 247)

(insn 25 24 26 2 ./CppUtilities.cc:247 (set (mem:QI (reg/v/f:DI 60 [ s ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 26 25 27 2 ("./CppUtilities.cc") 248)

(insn 27 26 28 2 ./CppUtilities.cc:248 (set (mem:QI (reg/v/f:DI 61 [ v ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(jump_insn 28 27 29 2 ./CppUtilities.cc:248 (set (pc)
        (label_ref:DI 127)) 548 {jump} (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 29 28 30)

;; Start of basic block 3, registers live: (nil)
(code_label 30 29 31 3 23 "" [1 uses])

(note 31 30 32 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 32 31 33 3 ("./CppUtilities.cc") 250)

(insn 33 32 34 3 ./CppUtilities.cc:250 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 62 [ signal ])
            (reg/v:DF 64 [ max ]))) -1 (nil)
    (nil))

(jump_insn 34 33 131 3 ./CppUtilities.cc:250 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 38)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 131 34 35 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(jump_insn 35 131 36 4 ./CppUtilities.cc:250 (set (pc)
        (label_ref 48)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

(barrier 36 35 38)

;; Start of basic block 5, registers live: (nil)
(code_label 38 36 39 5 29 "" [1 uses])

(note 39 38 40 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 40 39 41 5 ("./CppUtilities.cc") 252)

(insn 41 40 42 5 ./CppUtilities.cc:252 (set (mem:QI (reg/v/f:DI 59 [ h ]) [0 S1 A8])
        (const_int 120 [0x78])) -1 (nil)
    (nil))

(note 42 41 43 5 ("./CppUtilities.cc") 253)

(insn 43 42 44 5 ./CppUtilities.cc:253 (set (mem:QI (reg/v/f:DI 60 [ s ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 44 43 45 5 ("./CppUtilities.cc") 254)

(insn 45 44 46 5 ./CppUtilities.cc:254 (set (mem:QI (reg/v/f:DI 61 [ v ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(jump_insn 46 45 47 5 ./CppUtilities.cc:254 (set (pc)
        (label_ref:DI 127)) 548 {jump} (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 47 46 48)

;; Start of basic block 6, registers live: (nil)
(code_label 48 47 49 6 27 "" [1 uses])

(note 49 48 50 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 50 49 51 6 ("./CppUtilities.cc") 256)

(insn 51 50 52 6 ./CppUtilities.cc:256 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 62 [ signal ])
            (reg/v:DF 63 [ min ]))) -1 (nil)
    (nil))

(jump_insn 52 51 133 6 ./CppUtilities.cc:256 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 56)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 133 52 53 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(jump_insn 53 133 54 7 ./CppUtilities.cc:256 (set (pc)
        (label_ref 85)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

(barrier 54 53 56)

;; Start of basic block 8, registers live: (nil)
(code_label 56 54 57 8 32 "" [1 uses])

(note 57 56 58 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 58 57 59 8 ./CppUtilities.cc:256 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 65 [ med ])
            (reg/v:DF 62 [ signal ]))) -1 (nil)
    (nil))

(jump_insn 59 58 135 8 ./CppUtilities.cc:256 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 63)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 135 59 60 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(jump_insn 60 135 61 9 ./CppUtilities.cc:256 (set (pc)
        (label_ref 85)) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

(barrier 61 60 63)

;; Start of basic block 10, registers live: (nil)
(code_label 63 61 64 10 34 "" [1 uses])

(note 64 63 65 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 65 64 66 10 ("./CppUtilities.cc") 258)

(insn 66 65 67 10 ./CppUtilities.cc:258 (set (reg:DF 66)
        (minus:DF (reg/v:DF 65 [ med ])
            (reg/v:DF 63 [ min ]))) -1 (nil)
    (nil))

(insn 67 66 68 10 ./CppUtilities.cc:258 (set (reg:DF 68)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 68 67 69 10 ./CppUtilities.cc:258 (set (reg:DF 67)
        (div:DF (reg:DF 68)
            (reg:DF 66))) -1 (nil)
    (nil))

(insn 69 68 70 10 ./CppUtilities.cc:258 (set (reg:DF 69)
        (minus:DF (reg/v:DF 62 [ signal ])
            (reg/v:DF 63 [ min ]))) -1 (nil)
    (nil))

(insn 70 69 71 10 ./CppUtilities.cc:258 (set (reg/v:DF 58 [ t.267 ])
        (mult:DF (reg:DF 67)
            (reg:DF 69))) -1 (nil)
    (nil))

(note 71 70 72 10 ("./CppUtilities.cc") 260)

(insn 72 71 73 10 ./CppUtilities.cc:260 (set (reg:DF 70)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC4") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 73 72 74 10 ./CppUtilities.cc:260 (set (reg:DF 71)
        (mult:DF (reg/v:DF 58 [ t.267 ])
            (reg:DF 70))) -1 (nil)
    (nil))

(insn 74 73 75 10 ./CppUtilities.cc:260 (set (reg:DF 72)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC5") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.2e+2 [0x0.fp+7])
        (nil)))

(insn 75 74 76 10 ./CppUtilities.cc:260 (set (reg:DF 73)
        (mult:DF (reg:DF 71)
            (reg:DF 72))) -1 (nil)
    (nil))

(insn 76 75 77 10 ./CppUtilities.cc:260 (set (reg:SI 74)
        (fix:SI (reg:DF 73))) -1 (nil)
    (nil))

(insn 77 76 78 10 ./CppUtilities.cc:260 (set (mem:QI (reg/v/f:DI 59 [ h ]) [0 S1 A8])
        (subreg:QI (reg:SI 74) 0)) -1 (nil)
    (nil))

(note 78 77 79 10 ("./CppUtilities.cc") 261)

(insn 79 78 80 10 ./CppUtilities.cc:261 (set (mem:QI (reg/v/f:DI 60 [ s ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 80 79 81 10 ("./CppUtilities.cc") 262)

(insn 81 80 82 10 ./CppUtilities.cc:262 (set (mem:QI (reg/v/f:DI 61 [ v ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 82 81 83 10 ("./CppUtilities.cc") 256)

(jump_insn 83 82 84 10 ./CppUtilities.cc:256 (set (pc)
        (label_ref:DI 127)) 548 {jump} (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 84 83 85)

;; Start of basic block 11, registers live: (nil)
(code_label 85 84 86 11 30 "" [2 uses])

(note 86 85 87 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(note 87 86 88 11 ("./CppUtilities.cc") 264)

(insn 88 87 89 11 ./CppUtilities.cc:264 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 62 [ signal ])
            (reg/v:DF 65 [ med ]))) -1 (nil)
    (nil))

(jump_insn 89 88 137 11 ./CppUtilities.cc:264 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 93)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 6700 [0x1a2c])
        (nil)))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(note 137 89 90 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(jump_insn 90 137 91 12 ./CppUtilities.cc:264 (set (pc)
        (label_ref:DI 127)) 548 {jump} (nil)
    (nil))
;; End of basic block 12, registers live:
 (nil)

(barrier 91 90 93)

;; Start of basic block 13, registers live: (nil)
(code_label 93 91 94 13 36 "" [1 uses])

(note 94 93 95 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 95 94 96 13 ./CppUtilities.cc:264 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 64 [ max ])
            (reg/v:DF 62 [ signal ]))) -1 (nil)
    (nil))

(jump_insn 96 95 139 13 ./CppUtilities.cc:264 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 100)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 6700 [0x1a2c])
        (nil)))
;; End of basic block 13, registers live:
 (nil)

;; Start of basic block 14, registers live: (nil)
(note 139 96 97 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(jump_insn 97 139 98 14 ./CppUtilities.cc:264 (set (pc)
        (label_ref:DI 127)) 548 {jump} (nil)
    (nil))
;; End of basic block 14, registers live:
 (nil)

(barrier 98 97 100)

;; Start of basic block 15, registers live: (nil)
(code_label 100 98 101 15 38 "" [1 uses])

(note 101 100 102 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 102 101 103 15 ("./CppUtilities.cc") 268)

(insn 103 102 104 15 ./CppUtilities.cc:268 (set (reg:DF 75)
        (minus:DF (reg/v:DF 64 [ max ])
            (reg/v:DF 65 [ med ]))) -1 (nil)
    (nil))

(insn 104 103 105 15 ./CppUtilities.cc:268 (set (reg:DF 77)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 105 104 106 15 ./CppUtilities.cc:268 (set (reg:DF 76)
        (div:DF (reg:DF 77)
            (reg:DF 75))) -1 (nil)
    (nil))

(insn 106 105 107 15 ./CppUtilities.cc:268 (set (reg:DF 78)
        (minus:DF (reg/v:DF 62 [ signal ])
            (reg/v:DF 65 [ med ]))) -1 (nil)
    (nil))

(insn 107 106 108 15 ./CppUtilities.cc:268 (set (reg:DF 79)
        (mult:DF (reg:DF 76)
            (reg:DF 78))) -1 (nil)
    (nil))

(insn 108 107 109 15 ./CppUtilities.cc:268 (set (reg:DF 80)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC4") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 109 108 110 15 ./CppUtilities.cc:268 (set (reg:DF 81)
        (mult:DF (reg:DF 79)
            (reg:DF 80))) -1 (nil)
    (nil))

(insn 110 109 111 15 ./CppUtilities.cc:268 (set (reg:DF 82)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC4") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 111 110 112 15 ./CppUtilities.cc:268 (set (reg:DF 83)
        (plus:DF (reg:DF 81)
            (reg:DF 82))) -1 (nil)
    (nil))

(insn 112 111 113 15 ./CppUtilities.cc:268 (set (reg:DF 84)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC5") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.2e+2 [0x0.fp+7])
        (nil)))

(insn 113 112 114 15 ./CppUtilities.cc:268 (set (reg:DF 85)
        (mult:DF (reg:DF 83)
            (reg:DF 84))) -1 (nil)
    (nil))

(insn 114 113 115 15 ./CppUtilities.cc:268 (set (reg:SI 86)
        (fix:SI (reg:DF 85))) -1 (nil)
    (nil))

(insn 115 114 116 15 ./CppUtilities.cc:268 (set (mem:QI (reg/v/f:DI 59 [ h ]) [0 S1 A8])
        (subreg:QI (reg:SI 86) 0)) -1 (nil)
    (nil))

(note 116 115 117 15 ("./CppUtilities.cc") 269)

(insn 117 116 118 15 ./CppUtilities.cc:269 (set (mem:QI (reg/v/f:DI 60 [ s ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(note 118 117 119 15 ("./CppUtilities.cc") 270)

(insn 119 118 122 15 ./CppUtilities.cc:270 (set (mem:QI (reg/v/f:DI 61 [ v ]) [0 S1 A8])
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))
;; End of basic block 15, registers live:
 (nil)

(note 122 119 123 NOTE_INSN_FUNCTION_END)

(note 123 122 127 ("./CppUtilities.cc") 272)

;; Start of basic block 16, registers live: (nil)
(code_label 127 123 142 16 39 "" [5 uses])

(note 142 127 0 16 [bb 16] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 16, registers live:
 (nil)


;; Function void __static_initialization_and_destruction_0(int, int) (_Z41__static_initialization_and_destruction_0ii)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->4 redirected to 6
Deleted label in block 2.
Edge 2->4 redirected to 6
Deleted label in block 3.
Deleting block 4.
Deleting block 5.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("./CppUtilities.cc") 272)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppUtilities.cc:272 (set (reg/v:SI 58 [ __initialize_p ])
        (reg:SI 5 di [ __initialize_p ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:272 (set (reg/v:SI 59 [ __priority ])
        (reg:SI 4 si [ __priority ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppUtilities.cc") 272)

(insn 10 9 11 0 ./CppUtilities.cc:272 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 58 [ __initialize_p ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./CppUtilities.cc:272 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 33)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5467 [0x155b])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 13 11 14 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 14 13 15 1 ./CppUtilities.cc:272 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 59 [ __priority ])
            (const_int 65535 [0xffff]))) -1 (nil)
    (nil))

(jump_insn 15 14 17 1 ./CppUtilities.cc:272 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 33)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 7378 [0x1cd2])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 17 15 18 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 18 17 19 2 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

(insn 19 18 20 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2b4c434c5370 __ioinit>)) -1 (nil)
    (nil))

(call_insn 20 19 21 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitC1Ev") [flags 0x41] <function_decl 0x2b4c42c80b00 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 1 dx)
        (symbol_ref:DI ("__dso_handle") [flags 0x40] <var_decl 0x2b4c438edbb0 __dso_handle>)) -1 (nil)
    (nil))

(insn 22 21 23 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("__tcf_0") [flags 0x3] <function_decl 0x2b4c438ded00 __tcf_0>)) -1 (nil)
    (nil))

(call_insn/j 24 23 25 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_atexit") [flags 0x41] <function_decl 0x2b4c438dee00 __cxa_atexit>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 2, registers live:
 (nil)

(barrier 25 24 28)

(note 28 25 29 NOTE_INSN_FUNCTION_END)

(note 29 28 33 ("./CppUtilities.cc") 272)

;; Start of basic block 3, registers live: (nil)
(code_label 33 29 36 3 46 "" [2 uses])

(note 36 33 0 3 [bb 3] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 3, registers live:
 (nil)


;; Function (static initializers for ./CppUtilities.cc) (_GLOBAL__I__Z4joinPKcRSt6vectorIiSaIiEE)



try_optimize_cfg iteration 1

Deleting fallthru block 0.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 3 ("./CppUtilities.cc") 273)

(note 3 1 6 0 NOTE_INSN_FUNCTION_BEG)

;; Start of basic block 0, registers live: (nil)
(note 6 3 7 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(note 7 6 8 0 ("./CppUtilities.cc") 273)

(insn 8 7 9 0 ./CppUtilities.cc:273 (set (reg:SI 4 si)
        (const_int 65535 [0xffff])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppUtilities.cc:273 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 ./CppUtilities.cc:273 (call (mem:QI (symbol_ref:DI ("_Z41__static_initialization_and_destruction_0ii") [flags 0x3] <function_decl 0x2b4c438dec00 __static_initialization_and_destruction_0>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("./CppUtilities.cc") 273)


;; Function void __tcf_0(void*) (__tcf_0)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg/v/f:DI 58 [ D.37604 ])
        (reg:DI 5 di [ D.37604 ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

(insn 9 8 10 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2b4c434c5370 __ioinit>)) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitD1Ev") [flags 0x41] <function_decl 0x2b4c42c80d00 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)


;; Function void* SafeRealloc(void*, size_t) (_Z11SafeReallocPvm)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 3.
Merged 4 and 5 without moving.
Merged 4 and 6 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("./CppUtilities.cc") 118)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppUtilities.cc:118 (set (reg/v/f:DI 61 [ P ])
        (reg:DI 5 di [ P ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:118 (set (reg/v:DI 62 [ size ])
        (reg:DI 4 si [ size ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppUtilities.cc") 120)

(insn 10 9 11 0 ./CppUtilities.cc:120 (set (reg:DI 4 si)
        (reg/v:DI 62 [ size ])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppUtilities.cc:120 (set (reg:DI 5 di)
        (reg/v/f:DI 61 [ P ])) -1 (nil)
    (nil))

(call_insn 12 11 13 0 ./CppUtilities.cc:120 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("realloc") [flags 0x41] <function_decl 0x2b4c42429e00 realloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 13 12 14 0 ./CppUtilities.cc:120 (set (reg/f:DI 63)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 63)
        (nil)))

(insn 14 13 15 0 ./CppUtilities.cc:120 (set (reg/v/f:DI 58 [ Q ])
        (reg/f:DI 63)) -1 (nil)
    (nil))

(note 15 14 16 0 ("./CppUtilities.cc") 121)

(insn 16 15 17 0 ./CppUtilities.cc:121 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 58 [ Q ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 17 16 50 0 ./CppUtilities.cc:121 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 34)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8100 [0x1fa4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 50 17 18 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 18 50 19 1 ./CppUtilities.cc:121 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 62 [ size ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 19 18 21 1 ./CppUtilities.cc:121 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 34)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 21 19 22 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 22 21 23 2 ("./CppUtilities.cc") 122)

(insn 23 22 24 2 ./CppUtilities.cc:122 (set (reg:DI 5 di)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(call_insn 24 23 25 2 ./CppUtilities.cc:122 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_allocate_exception") [flags 0x41] <function_decl 0x2b4c43769c00 __cxa_allocate_exception>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 25 24 26 2 ./CppUtilities.cc:122 (set (reg:DI 59 [ D.36102 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 26 25 27 2 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new") 58)

(insn 27 26 28 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new:58 (set (mem/s:DI (reg:DI 59 [ D.36102 ]) [6 <variable>.D.10655._vptr.exception+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259dd10 _ZTVSt9bad_alloc>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 28 27 29 2 ("./CppUtilities.cc") 122)

(insn 29 28 30 2 ./CppUtilities.cc:122 (set (reg:DI 1 dx)
        (symbol_ref:DI ("_ZNSt9bad_allocD1Ev") [flags 0x41] <function_decl 0x2b4c4259f200 __comp_dtor >)) -1 (nil)
    (nil))

(insn 30 29 31 2 ./CppUtilities.cc:122 (set (reg:DI 4 si)
        (symbol_ref:DI ("_ZTISt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259ddc0 _ZTISt9bad_alloc>)) -1 (nil)
    (nil))

(insn 31 30 32 2 ./CppUtilities.cc:122 (set (reg:DI 5 di)
        (reg:DI 59 [ D.36102 ])) -1 (nil)
    (nil))

(call_insn 32 31 33 2 ./CppUtilities.cc:122 (call (mem:QI (symbol_ref:DI ("__cxa_throw") [flags 0x41] <function_decl 0x2b4c43769b00 __cxa_throw>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 2, registers live:
 (nil)

(barrier 33 32 34)

;; Start of basic block 3, registers live: (nil)
(code_label 34 33 35 3 53 "" [2 uses])

(note 35 34 36 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 36 35 39 3 ./CppUtilities.cc:122 (set (reg:DI 60 [ <result> ])
        (reg/v/f:DI 58 [ Q ])) -1 (nil)
    (nil))

(note 39 36 40 3 NOTE_INSN_FUNCTION_END)

(note 40 39 42 3 ("./CppUtilities.cc") 124)

(insn 42 40 48 3 ./CppUtilities.cc:124 (set (reg/i:DI 0 ax)
        (reg:DI 60 [ <result> ])) -1 (nil)
    (nil))

(insn 48 42 0 3 ./CppUtilities.cc:124 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)


;; Function void* SafeMalloc(size_t) (_Z10SafeMallocm)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 3 and 4 without moving.
Merged 3 and 5 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./CppUtilities.cc") 108)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./CppUtilities.cc:108 (set (reg/v:DI 61 [ size ])
        (reg:DI 5 di [ size ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./CppUtilities.cc") 110)

(insn 9 8 10 0 ./CppUtilities.cc:110 (set (reg:DI 5 di)
        (reg/v:DI 61 [ size ])) -1 (nil)
    (nil))

(call_insn 10 9 11 0 ./CppUtilities.cc:110 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b4c40f29100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 11 10 12 0 ./CppUtilities.cc:110 (set (reg/f:DI 62)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 62)
        (nil)))

(insn 12 11 13 0 ./CppUtilities.cc:110 (set (reg/v/f:DI 58 [ Q ])
        (reg/f:DI 62)) -1 (nil)
    (nil))

(note 13 12 14 0 ("./CppUtilities.cc") 111)

(insn 14 13 15 0 ./CppUtilities.cc:111 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 58 [ Q ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 15 14 17 0 ./CppUtilities.cc:111 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 30)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 17 15 18 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 18 17 19 1 ("./CppUtilities.cc") 112)

(insn 19 18 20 1 ./CppUtilities.cc:112 (set (reg:DI 5 di)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(call_insn 20 19 21 1 ./CppUtilities.cc:112 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_allocate_exception") [flags 0x41] <function_decl 0x2b4c43769c00 __cxa_allocate_exception>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 1 ./CppUtilities.cc:112 (set (reg:DI 59 [ D.36091 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 22 21 23 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new") 58)

(insn 23 22 24 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new:58 (set (mem/s:DI (reg:DI 59 [ D.36091 ]) [6 <variable>.D.10655._vptr.exception+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259dd10 _ZTVSt9bad_alloc>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 24 23 25 1 ("./CppUtilities.cc") 112)

(insn 25 24 26 1 ./CppUtilities.cc:112 (set (reg:DI 1 dx)
        (symbol_ref:DI ("_ZNSt9bad_allocD1Ev") [flags 0x41] <function_decl 0x2b4c4259f200 __comp_dtor >)) -1 (nil)
    (nil))

(insn 26 25 27 1 ./CppUtilities.cc:112 (set (reg:DI 4 si)
        (symbol_ref:DI ("_ZTISt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259ddc0 _ZTISt9bad_alloc>)) -1 (nil)
    (nil))

(insn 27 26 28 1 ./CppUtilities.cc:112 (set (reg:DI 5 di)
        (reg:DI 59 [ D.36091 ])) -1 (nil)
    (nil))

(call_insn 28 27 29 1 ./CppUtilities.cc:112 (call (mem:QI (symbol_ref:DI ("__cxa_throw") [flags 0x41] <function_decl 0x2b4c43769b00 __cxa_throw>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 1, registers live:
 (nil)

(barrier 29 28 30)

;; Start of basic block 2, registers live: (nil)
(code_label 30 29 31 2 58 "" [1 uses])

(note 31 30 32 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 32 31 35 2 ./CppUtilities.cc:112 (set (reg:DI 60 [ <result> ])
        (reg/v/f:DI 58 [ Q ])) -1 (nil)
    (nil))

(note 35 32 36 2 NOTE_INSN_FUNCTION_END)

(note 36 35 38 2 ("./CppUtilities.cc") 114)

(insn 38 36 44 2 ./CppUtilities.cc:114 (set (reg/i:DI 0 ax)
        (reg:DI 60 [ <result> ])) -1 (nil)
    (nil))

(insn 44 38 0 2 ./CppUtilities.cc:114 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)


;; Function bool ReadLine(FILE*, char*&, size_t&) (_Z8ReadLineP8_IO_FILERPcRm)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 4 failed.
Forwarding edge 0->2 to 4 failed.
Deleted label in block 5.
Forwarding edge 5->6 to 8 failed.
Merged 9 and 10 without moving.
Merged 9 and 11 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 4 failed.
Forwarding edge 5->6 to 8 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 3 failed.
Forwarding edge 4->5 to 7 failed.
(note 1 0 8 ("./CppUtilities.cc") 54)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./CppUtilities.cc:54 (set (reg/v/f:DI 69 [ file ])
        (reg:DI 5 di [ file ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:54 (set (reg/v/f:DI 70 [ buff ])
        (reg:DI 4 si [ buff ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppUtilities.cc:54 (set (reg/v/f:DI 71 [ size ])
        (reg:DI 1 dx [ size ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./CppUtilities.cc") 56)

(insn 11 10 12 0 ./CppUtilities.cc:56 (set (reg:DI 72)
        (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppUtilities.cc:56 (set (reg:DI 73)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppUtilities.cc:56 (set (mem:QI (plus:DI (plus:DI (reg:DI 72)
                    (reg:DI 73))
                (const_int -1 [0xffffffffffffffff])) [0 S1 A8])
        (const_int 36 [0x24])) -1 (nil)
    (nil))

(note 14 13 15 0 ("./CppUtilities.cc") 57)

(insn 15 14 16 0 ./CppUtilities.cc:57 (set (reg:DI 74)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppUtilities.cc:57 (set (reg:DI 1 dx)
        (reg/v/f:DI 69 [ file ])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppUtilities.cc:57 (set (reg:SI 4 si)
        (subreg:SI (reg:DI 74) 0)) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppUtilities.cc:57 (set (reg:DI 5 di)
        (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])) -1 (nil)
    (nil))

(call_insn 19 18 20 0 ./CppUtilities.cc:57 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fgets") [flags 0x41] <function_decl 0x2b4c4221a800 fgets>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 20 19 21 0 ./CppUtilities.cc:57 (set (reg:DI 67 [ D.35969 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppUtilities.cc:57 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 67 [ D.35969 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 22 21 25 0 ./CppUtilities.cc:57 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1036 [0x40c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 25 22 23 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 23 25 24 1 ./CppUtilities.cc:57 (set (pc)
        (label_ref 56)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 24 23 26)

;; Start of basic block 2, registers live: (nil)
(code_label 26 24 27 2 65 "" [1 uses])

(note 27 26 28 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 28 27 29 2 ("./CppUtilities.cc") 62)

(insn 29 28 30 2 ./CppUtilities.cc:62 (set (reg/v:DI 63 [ old ])
        (reg:DI 60 [ temp.442 ])) -1 (nil)
    (nil))

(note 30 29 31 2 ("./CppUtilities.cc") 63)

(insn 31 30 32 2 ./CppUtilities.cc:63 (parallel [
            (set (reg:DI 59 [ temp.443 ])
                (ashift:DI (reg:DI 60 [ temp.442 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 32 31 33 2 ./CppUtilities.cc:63 (set (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])
        (reg:DI 59 [ temp.443 ])) -1 (nil)
    (nil))

(note 33 32 34 2 ("./CppUtilities.cc") 64)

(insn 34 33 35 2 ./CppUtilities.cc:64 (set (reg:DI 4 si)
        (reg:DI 59 [ temp.443 ])) -1 (nil)
    (nil))

(insn 35 34 36 2 ./CppUtilities.cc:64 (set (reg:DI 5 di)
        (reg:DI 62 [ temp.438 ])) -1 (nil)
    (nil))

(call_insn 36 35 37 2 ./CppUtilities.cc:64 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Z11SafeReallocPvm") [flags 0x3] <function_decl 0x2b4c4365c300 SafeRealloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 37 36 38 2 ./CppUtilities.cc:64 (set (reg:DI 65 [ D.35979 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 38 37 39 2 ./CppUtilities.cc:64 (set (reg:DI 61 [ temp.439 ])
        (reg:DI 65 [ D.35979 ])) -1 (nil)
    (nil))

(insn 39 38 40 2 ./CppUtilities.cc:64 (set (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])
        (reg:DI 61 [ temp.439 ])) -1 (nil)
    (nil))

(note 40 39 41 2 ("./CppUtilities.cc") 66)

(insn 41 40 42 2 ./CppUtilities.cc:66 (set (reg:DI 75)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 42 41 43 2 ./CppUtilities.cc:66 (set (mem:QI (plus:DI (plus:DI (reg:DI 61 [ temp.439 ])
                    (reg:DI 75))
                (const_int -1 [0xffffffffffffffff])) [0 S1 A8])
        (const_int 36 [0x24])) -1 (nil)
    (nil))

(note 43 42 44 2 ("./CppUtilities.cc") 67)

(insn 44 43 45 2 ./CppUtilities.cc:67 (set (reg:DI 76)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 45 44 46 2 ./CppUtilities.cc:67 (parallel [
            (set (reg:SI 77)
                (plus:SI (subreg:SI (reg:DI 76) 0)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 46 45 47 2 ./CppUtilities.cc:67 (parallel [
            (set (reg:SI 78)
                (minus:SI (reg:SI 77)
                    (subreg:SI (reg/v:DI 63 [ old ]) 0)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 47 46 48 2 ./CppUtilities.cc:67 (parallel [
            (set (reg:DI 79)
                (plus:DI (reg/v:DI 63 [ old ])
                    (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 48 47 49 2 ./CppUtilities.cc:67 (parallel [
            (set (reg:DI 80)
                (plus:DI (reg:DI 79)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 49 48 50 2 ./CppUtilities.cc:67 (set (reg:DI 1 dx)
        (reg/v/f:DI 69 [ file ])) -1 (nil)
    (nil))

(insn 50 49 51 2 ./CppUtilities.cc:67 (set (reg:SI 4 si)
        (reg:SI 78)) -1 (nil)
    (nil))

(insn 51 50 52 2 ./CppUtilities.cc:67 (set (reg:DI 5 di)
        (reg:DI 80)) -1 (nil)
    (nil))

(call_insn 52 51 53 2 ./CppUtilities.cc:67 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fgets") [flags 0x41] <function_decl 0x2b4c4221a800 fgets>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 53 52 54 2 ./CppUtilities.cc:67 (set (reg:DI 64 [ D.35991 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 54 53 55 2 ./CppUtilities.cc:67 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 64 [ D.35991 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 55 54 56 2 ./CppUtilities.cc:67 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1100 [0x44c])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 56 55 57 3 64 "" [1 uses])

(note 57 56 58 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 58 57 59 3 ("./CppUtilities.cc") 60)

(insn 59 58 60 3 ./CppUtilities.cc:60 (set (reg:DI 62 [ temp.438 ])
        (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])) -1 (nil)
    (nil))

(insn 60 59 61 3 ./CppUtilities.cc:60 (set (reg:DI 60 [ temp.442 ])
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 61 60 62 3 ./CppUtilities.cc:60 (set (reg:DI 58 [ temp.448 ])
        (reg:DI 60 [ temp.442 ])) -1 (nil)
    (nil))

(insn 62 61 63 3 ./CppUtilities.cc:60 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:QI (plus:DI (plus:DI (reg:DI 62 [ temp.438 ])
                        (reg:DI 58 [ temp.448 ]))
                    (const_int -1 [0xffffffffffffffff])) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 63 62 65 3 ./CppUtilities.cc:60 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 76)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 65 63 66 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 66 65 67 4 ./CppUtilities.cc:60 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:QI (plus:DI (plus:DI (reg:DI 62 [ temp.438 ])
                        (reg:DI 58 [ temp.448 ]))
                    (const_int -2 [0xfffffffffffffffe])) [0 S1 A8])
            (const_int 10 [0xa]))) -1 (nil)
    (nil))

(jump_insn 67 66 70 4 ./CppUtilities.cc:60 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 26)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 70 67 68 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 68 70 69 5 ./CppUtilities.cc:60 (set (pc)
        (label_ref 76)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 69 68 71)

;; Start of basic block 6, registers live: (nil)
(code_label 71 69 72 6 63 "" [2 uses])

(note 72 71 73 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 73 72 74 6 ./CppUtilities.cc:60 (set (reg:SI 66 [ D.35970 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 74 73 75 6 ./CppUtilities.cc:60 (set (pc)
        (label_ref 79)) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

(barrier 75 74 76)

;; Start of basic block 7, registers live: (nil)
(code_label 76 75 77 7 66 "" [2 uses])

(note 77 76 78 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 78 77 79 7 ./CppUtilities.cc:60 (set (reg:SI 66 [ D.35970 ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 79 78 80 8 68 "" [1 uses])

(note 80 79 81 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 81 80 84 8 ./CppUtilities.cc:60 (set (reg:SI 68 [ <result> ])
        (reg:SI 66 [ D.35970 ])) -1 (nil)
    (nil))

(note 84 81 85 8 NOTE_INSN_FUNCTION_END)

(note 85 84 87 8 ("./CppUtilities.cc") 72)

(insn 87 85 93 8 ./CppUtilities.cc:72 (set (reg/i:SI 0 ax)
        (reg:SI 68 [ <result> ])) -1 (nil)
    (nil))

(insn 93 87 0 8 ./CppUtilities.cc:72 (use (reg/i:SI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)


;; Function static _CharT* std::basic_string<_CharT, _Traits, _Alloc>::_S_construct(_InIterator, _InIterator, const _Alloc&, std::forward_iterator_tag) [with _FwdIterator = char*, _CharT = char, _Traits = std::char_traits<char>, _Alloc = std::allocator<char>] (_ZNSs12_S_constructIPcEES0_T_S1_RKSaIcESt20forward_iterator_tag)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 5.
Forwarding edge 9->10 to 4 failed.
Merged 11 and 12 without moving.


try_optimize_cfg iteration 2

Forwarding edge 9->10 to 4 failed.


try_optimize_cfg iteration 1

Forwarding edge 8->9 to 3 failed.
(note 1 0 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 137)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:137 (set (reg/v/f:DI 66 [ __beg ])
        (reg:DI 5 di [ __beg ])) -1 (nil)
    (nil))

(insn 4 3 5 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:137 (set (reg/v/f:DI 67 [ __end ])
        (reg:DI 4 si [ __end ])) -1 (nil)
    (nil))

(insn 5 4 6 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:137 (set (reg/v/f:DI 68 [ __a ])
        (reg:DI 1 dx [ __a ])) -1 (nil)
    (nil))

(insn 6 5 7 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:137 (set (reg/v:QI 69 [ D.37324 ])
        (mem/s/c:QI (reg/f:DI 53 virtual-incoming-args) [28 D.37324+0 S1 A64])) -1 (nil)
    (expr_list:REG_EQUIV (mem/s/c:QI (reg/f:DI 53 virtual-incoming-args) [28 D.37324+0 S1 A64])
        (nil)))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 140)

(insn 12 11 13 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:140 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 66 [ __beg ])
            (reg/v/f:DI 67 [ __end ]))) -1 (nil)
    (nil))

(jump_insn 13 12 15 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:140 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 81)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 6463 [0x193f])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 15 13 16 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 16 15 17 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 17 16 18 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 62 [ this ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 18 17 19 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 213)

(insn 19 18 20 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:213 (parallel [
            (set (reg/v/f:DI 60 [ __s1 ])
                (plus:DI (reg:DI 62 [ this ])
                    (const_int 24 [0x18])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 20 19 21 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 141)

(jump_insn 21 20 22 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:141 (set (pc)
        (label_ref 76)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 22 21 23)

;; Start of basic block 2, registers live: (nil)
(code_label 23 22 24 2 75 "" [1 uses])

(note 24 23 25 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 25 24 26 2 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 145)

(insn 26 25 27 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:145 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC6") [flags 0x2] <string_cst 0x2b4c43ae1c40>)) -1 (nil)
    (nil))

(call_insn 27 26 28 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:145 (call (mem:QI (symbol_ref:DI ("_ZSt19__throw_logic_errorPKc") [flags 0x41] <function_decl 0x2b4c423c1c00 __throw_logic_error>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 2, registers live:
 (nil)

(barrier 28 27 29)

;; Start of basic block 3, registers live: (nil)
(code_label 29 28 30 3 76 "" [1 uses])

(note 30 29 31 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 31 30 32 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator_base_funcs.h") 97)

(insn 32 31 33 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator_base_funcs.h:97 (parallel [
            (set (reg:DI 61 [ D.38268 ])
                (minus:DI (reg/v/f:DI 67 [ __end ])
                    (reg/v/f:DI 66 [ __beg ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 33 32 34 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 148)

(insn 34 33 35 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:148 (set (reg/v:DI 64 [ __dnew ])
        (reg:DI 61 [ D.38268 ])) -1 (nil)
    (nil))

(note 35 34 36 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 150)

(insn 36 35 37 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:150 (set (reg:DI 1 dx)
        (reg/v/f:DI 68 [ __a ])) -1 (nil)
    (nil))

(insn 37 36 38 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:150 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 38 37 39 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:150 (set (reg:DI 5 di)
        (reg/v:DI 64 [ __dnew ])) -1 (nil)
    (nil))

(call_insn 39 38 40 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:150 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep9_S_createEmmRKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2e900 _S_create>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 40 39 41 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:150 (set (reg/v/f:DI 63 [ __r ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 41 40 42 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 213)

(insn 42 41 43 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:213 (parallel [
            (set (reg/v/f:DI 60 [ __s1 ])
                (plus:DI (reg/v/f:DI 63 [ __r ])
                    (const_int 24 [0x18])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 43 42 44 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:213 (set (reg/v/f:DI 59 [ __s2 ])
        (reg/v/f:DI 66 [ __beg ])) -1 (nil)
    (nil))

(note 44 43 45 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 338)

(insn 45 44 46 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:338 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 61 [ D.38268 ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 46 45 48 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:338 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 54)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5120 [0x1400])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 48 46 49 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 49 48 50 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h") 241)

(insn 50 49 51 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:241 (set (reg:QI 70)
        (mem:QI (reg/v/f:DI 59 [ __s2 ]) [0 S1 A8])) -1 (nil)
    (nil))

(insn 51 50 52 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:241 (set (mem:QI (reg/v/f:DI 60 [ __s1 ]) [0 S1 A8])
        (reg:QI 70)) -1 (nil)
    (nil))

(jump_insn 52 51 53 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:241 (set (pc)
        (label_ref 68)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

(barrier 53 52 54)

;; Start of basic block 5, registers live: (nil)
(code_label 54 53 55 5 77 "" [1 uses])

(note 55 54 56 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 56 55 59 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h") 269)

(insn 59 56 60 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 73)
        (reg/v/f:DI 60 [ __s1 ])) -1 (nil)
    (nil))

(insn 60 59 61 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 74)
        (reg/v/f:DI 59 [ __s2 ])) -1 (nil)
    (nil))

(insn 61 60 62 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 75)
        (reg/v:DI 64 [ __dnew ])) -1 (nil)
    (nil))

(insn 62 61 63 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 1 dx)
        (reg:DI 75)) -1 (nil)
    (nil))

(insn 63 62 64 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 4 si)
        (reg:DI 74)) -1 (nil)
    (nil))

(insn 64 63 65 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 5 di)
        (reg:DI 73)) -1 (nil)
    (nil))

(call_insn 65 64 66 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memcpy") [flags 0x41] <function_decl 0x2b4c43aec600 memcpy>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 66 65 67 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 76)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 67 66 68 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:269 (set (reg:DI 58 [ D.38957 ])
        (reg:DI 76)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(code_label 68 67 69 6 79 "" [1 uses])

(note 69 68 70 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 70 69 71 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 199)

(insn 71 70 72 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:199 (set (mem/s:SI (plus:DI (reg/v/f:DI 63 [ __r ])
                (const_int 16 [0x10])) [16 <variable>.D.17834._M_refcount+0 S4 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 72 71 73 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 205)

(insn 73 72 74 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:205 (set (mem/s:DI (reg/v/f:DI 63 [ __r ]) [19 <variable>.D.17834._M_length+0 S8 A64])
        (reg/v:DI 64 [ __dnew ])) -1 (nil)
    (nil))

(note 74 73 75 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h") 241)

(insn 75 74 76 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:241 (set (mem:QI (plus:DI (reg/v/f:DI 60 [ __s1 ])
                (reg/v:DI 64 [ __dnew ])) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(code_label 76 75 77 7 74 "" [1 uses])

(note 77 76 78 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 78 77 79 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:241 (set (reg:DI 65 [ <result> ])
        (reg/v/f:DI 60 [ __s1 ])) -1 (nil)
    (nil))

(jump_insn 79 78 80 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/char_traits.h:241 (set (pc)
        (label_ref 91)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

(barrier 80 79 81)

;; Start of basic block 8, registers live: (nil)
(code_label 81 80 82 8 72 "" [1 uses])

(note 82 81 83 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 83 82 84 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 144)

(insn 84 83 85 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:144 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 66 [ __beg ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 85 84 88 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:144 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 23)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 99 [0x63])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 99 [0x63]))
            (nil))))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 88 85 86 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(jump_insn 86 88 87 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:144 (set (pc)
        (label_ref 29)) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

(barrier 87 86 89)

(note 89 87 90 NOTE_INSN_FUNCTION_END)

(note 90 89 91 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 159)

;; Start of basic block 10, registers live: (nil)
(code_label 91 90 100 10 71 "" [1 uses])

(note 100 91 92 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 92 100 98 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:159 (set (reg/i:DI 0 ax)
        (reg:DI 65 [ <result> ])) -1 (nil)
    (nil))

(insn 98 92 0 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:159 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)


;; Function std::basic_string<_CharT, _Traits, _Alloc>::basic_string(_InputIterator, _InputIterator, const _Alloc&) [with _InputIterator = char*, _CharT = char, _Traits = std::char_traits<char>, _Alloc = std::allocator<char>] (_ZNSsC1IPcEET_S1_RKSaIcE)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 3 failed.
Merged 0 and 2 without moving.
Merged 0 and 3 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 9 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 240)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:240 (set (reg/f:DI 60 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(insn 4 3 5 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:240 (set (reg/v/f:DI 61 [ __beg ])
        (reg:DI 4 si [ __beg ])) -1 (nil)
    (nil))

(insn 5 4 6 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:240 (set (reg/v/f:DI 62 [ __end ])
        (reg:DI 1 dx [ __end ])) -1 (nil)
    (nil))

(insn 6 5 7 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc:240 (set (reg/v/f:DI 63 [ __a ])
        (reg:DI 2 cx [ __a ])) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1449)

(insn 12 11 13 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1449 (set (mem:QI (reg/f:DI 56 virtual-outgoing-args) [0 S1 A64])
        (reg:QI 58 [ D.39025 ])) -1 (nil)
    (nil))

(insn 13 12 14 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1449 (set (reg:DI 1 dx)
        (reg/v/f:DI 63 [ __a ])) -1 (nil)
    (nil))

(insn 14 13 15 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1449 (set (reg:DI 4 si)
        (reg/v/f:DI 62 [ __end ])) -1 (nil)
    (nil))

(insn 15 14 16 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1449 (set (reg:DI 5 di)
        (reg/v/f:DI 61 [ __beg ])) -1 (nil)
    (nil))

(call_insn 16 15 17 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1449 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref/i:DI ("_ZNSs12_S_constructIPcEES0_T_S1_RKSaIcESt20forward_iterator_tag") [flags 0x1] <function_decl 0x2b4c4389ca00 _S_construct>) [0 S1 A8])
            (const_int 8 [0x8]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 17 16 18 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1449 (set (reg/v/f:DI 59 [ __dat ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 18 17 19 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 257)

(insn 19 18 20 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:257 (set (mem/s:DI (reg/f:DI 60 [ this ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])
        (reg/v/f:DI 59 [ __dat ])) -1 (nil)
    (nil))
;; End of basic block 0, registers live:
 (nil)

(note 20 19 21 NOTE_INSN_FUNCTION_END)

(note 21 20 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.tcc") 242)


;; Function void* SafeCalloc(size_t, size_t) (_Z10SafeCallocmm)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 3 and 4 without moving.
Merged 3 and 5 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("./CppUtilities.cc") 98)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppUtilities.cc:98 (set (reg/v:DI 61 [ num ])
        (reg:DI 5 di [ num ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:98 (set (reg/v:DI 62 [ size ])
        (reg:DI 4 si [ size ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppUtilities.cc") 100)

(insn 10 9 11 0 ./CppUtilities.cc:100 (set (reg:DI 4 si)
        (reg/v:DI 62 [ size ])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppUtilities.cc:100 (set (reg:DI 5 di)
        (reg/v:DI 61 [ num ])) -1 (nil)
    (nil))

(call_insn 12 11 13 0 ./CppUtilities.cc:100 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b4c40f18200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 13 12 14 0 ./CppUtilities.cc:100 (set (reg/f:DI 63)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 63)
        (nil)))

(insn 14 13 15 0 ./CppUtilities.cc:100 (set (reg/v/f:DI 58 [ Q ])
        (reg/f:DI 63)) -1 (nil)
    (nil))

(note 15 14 16 0 ("./CppUtilities.cc") 101)

(insn 16 15 17 0 ./CppUtilities.cc:101 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 58 [ Q ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 17 16 19 0 ./CppUtilities.cc:101 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 32)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 19 17 20 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 20 19 21 1 ("./CppUtilities.cc") 102)

(insn 21 20 22 1 ./CppUtilities.cc:102 (set (reg:DI 5 di)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(call_insn 22 21 23 1 ./CppUtilities.cc:102 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_allocate_exception") [flags 0x41] <function_decl 0x2b4c43769c00 __cxa_allocate_exception>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 23 22 24 1 ./CppUtilities.cc:102 (set (reg:DI 59 [ D.36059 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 24 23 25 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new") 58)

(insn 25 24 26 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new:58 (set (mem/s:DI (reg:DI 59 [ D.36059 ]) [6 <variable>.D.10655._vptr.exception+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259dd10 _ZTVSt9bad_alloc>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 26 25 27 1 ("./CppUtilities.cc") 102)

(insn 27 26 28 1 ./CppUtilities.cc:102 (set (reg:DI 1 dx)
        (symbol_ref:DI ("_ZNSt9bad_allocD1Ev") [flags 0x41] <function_decl 0x2b4c4259f200 __comp_dtor >)) -1 (nil)
    (nil))

(insn 28 27 29 1 ./CppUtilities.cc:102 (set (reg:DI 4 si)
        (symbol_ref:DI ("_ZTISt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259ddc0 _ZTISt9bad_alloc>)) -1 (nil)
    (nil))

(insn 29 28 30 1 ./CppUtilities.cc:102 (set (reg:DI 5 di)
        (reg:DI 59 [ D.36059 ])) -1 (nil)
    (nil))

(call_insn 30 29 31 1 ./CppUtilities.cc:102 (call (mem:QI (symbol_ref:DI ("__cxa_throw") [flags 0x41] <function_decl 0x2b4c43769b00 __cxa_throw>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 1, registers live:
 (nil)

(barrier 31 30 32)

;; Start of basic block 2, registers live: (nil)
(code_label 32 31 33 2 87 "" [1 uses])

(note 33 32 34 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 34 33 37 2 ./CppUtilities.cc:102 (set (reg:DI 60 [ <result> ])
        (reg/v/f:DI 58 [ Q ])) -1 (nil)
    (nil))

(note 37 34 38 2 NOTE_INSN_FUNCTION_END)

(note 38 37 40 2 ("./CppUtilities.cc") 104)

(insn 40 38 46 2 ./CppUtilities.cc:104 (set (reg/i:DI 0 ax)
        (reg:DI 60 [ <result> ])) -1 (nil)
    (nil))

(insn 46 40 0 2 ./CppUtilities.cc:104 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)


;; Function char* SafeStrdup(const char*) (_Z10SafeStrdupPKc)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 3 and 4 without moving.
Merged 3 and 5 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./CppUtilities.cc") 128)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./CppUtilities.cc:128 (set (reg/v/f:DI 61 [ str ])
        (reg:DI 5 di [ str ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./CppUtilities.cc") 130)

(insn 9 8 10 0 ./CppUtilities.cc:130 (set (reg:DI 5 di)
        (reg/v/f:DI 61 [ str ])) -1 (nil)
    (nil))

(call_insn 10 9 11 0 ./CppUtilities.cc:130 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("strdup") [flags 0x41] <function_decl 0x2b4c40ef4c00 strdup>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 11 10 12 0 ./CppUtilities.cc:130 (set (reg/f:DI 62)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 62)
        (nil)))

(insn 12 11 13 0 ./CppUtilities.cc:130 (set (reg/v/f:DI 58 [ Q ])
        (reg/f:DI 62)) -1 (nil)
    (nil))

(note 13 12 14 0 ("./CppUtilities.cc") 131)

(insn 14 13 15 0 ./CppUtilities.cc:131 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v/f:DI 58 [ Q ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 15 14 17 0 ./CppUtilities.cc:131 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 30)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 17 15 18 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 18 17 19 1 ("./CppUtilities.cc") 132)

(insn 19 18 20 1 ./CppUtilities.cc:132 (set (reg:DI 5 di)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(call_insn 20 19 21 1 ./CppUtilities.cc:132 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_allocate_exception") [flags 0x41] <function_decl 0x2b4c43769c00 __cxa_allocate_exception>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 1 ./CppUtilities.cc:132 (set (reg:DI 59 [ D.36115 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 22 21 23 1 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new") 58)

(insn 23 22 24 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/new:58 (set (mem/s:DI (reg:DI 59 [ D.36115 ]) [6 <variable>.D.10655._vptr.exception+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259dd10 _ZTVSt9bad_alloc>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 24 23 25 1 ("./CppUtilities.cc") 132)

(insn 25 24 26 1 ./CppUtilities.cc:132 (set (reg:DI 1 dx)
        (symbol_ref:DI ("_ZNSt9bad_allocD1Ev") [flags 0x41] <function_decl 0x2b4c4259f200 __comp_dtor >)) -1 (nil)
    (nil))

(insn 26 25 27 1 ./CppUtilities.cc:132 (set (reg:DI 4 si)
        (symbol_ref:DI ("_ZTISt9bad_alloc") [flags 0x40] <var_decl 0x2b4c4259ddc0 _ZTISt9bad_alloc>)) -1 (nil)
    (nil))

(insn 27 26 28 1 ./CppUtilities.cc:132 (set (reg:DI 5 di)
        (reg:DI 59 [ D.36115 ])) -1 (nil)
    (nil))

(call_insn 28 27 29 1 ./CppUtilities.cc:132 (call (mem:QI (symbol_ref:DI ("__cxa_throw") [flags 0x41] <function_decl 0x2b4c43769b00 __cxa_throw>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 1, registers live:
 (nil)

(barrier 29 28 30)

;; Start of basic block 2, registers live: (nil)
(code_label 30 29 31 2 92 "" [1 uses])

(note 31 30 32 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 32 31 35 2 ./CppUtilities.cc:132 (set (reg:DI 60 [ <result> ])
        (reg/v/f:DI 58 [ Q ])) -1 (nil)
    (nil))

(note 35 32 36 2 NOTE_INSN_FUNCTION_END)

(note 36 35 38 2 ("./CppUtilities.cc") 134)

(insn 38 36 44 2 ./CppUtilities.cc:134 (set (reg/i:DI 0 ax)
        (reg:DI 60 [ <result> ])) -1 (nil)
    (nil))

(insn 44 38 0 2 ./CppUtilities.cc:134 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)


;; Function bool ReadLine(std::istream&, char*&, size_t&) (_Z8ReadLineRSiRPcRm)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 4 failed.
Forwarding edge 0->2 to 4 failed.
Deleted label in block 5.
Forwarding edge 5->6 to 8 failed.
Merged 9 and 10 without moving.
Merged 9 and 11 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 4 failed.
Forwarding edge 5->6 to 8 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 3 failed.
Forwarding edge 4->5 to 7 failed.
(note 1 0 8 ("./CppUtilities.cc") 76)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./CppUtilities.cc:76 (set (reg/v/f:DI 69 [ file ])
        (reg:DI 5 di [ file ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:76 (set (reg/v/f:DI 70 [ buff ])
        (reg:DI 4 si [ buff ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppUtilities.cc:76 (set (reg/v/f:DI 71 [ size ])
        (reg:DI 1 dx [ size ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./CppUtilities.cc") 78)

(insn 11 10 12 0 ./CppUtilities.cc:78 (set (reg:DI 72)
        (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppUtilities.cc:78 (set (reg:DI 73)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppUtilities.cc:78 (set (mem:QI (plus:DI (plus:DI (reg:DI 72)
                    (reg:DI 73))
                (const_int -1 [0xffffffffffffffff])) [0 S1 A8])
        (const_int 36 [0x24])) -1 (nil)
    (nil))

(note 14 13 15 0 ("./CppUtilities.cc") 79)

(insn 15 14 16 0 ./CppUtilities.cc:79 (set (reg:DI 1 dx)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppUtilities.cc:79 (set (reg:DI 4 si)
        (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppUtilities.cc:79 (set (reg:DI 5 di)
        (reg/v/f:DI 69 [ file ])) -1 (nil)
    (nil))

(call_insn 18 17 19 0 ./CppUtilities.cc:79 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSi4readEPcl") [flags 0x41] <function_decl 0x2b4c433c7400 read>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 19 18 20 0 ./CppUtilities.cc:79 (set (reg/f:DI 67 [ D.36005 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 20 19 21 0 ./CppUtilities.cc:79 (set (reg/f:DI 74)
        (mem/s/f:DI (reg/f:DI 67 [ D.36005 ]) [6 <variable>._vptr.basic_istream+0 S8 A64])) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppUtilities.cc:79 (set (reg:DI 75)
        (mem:DI (plus:DI (reg/f:DI 74)
                (const_int -24 [0xffffffffffffffe8])) [19 S8 A64])) -1 (nil)
    (nil))

(insn 22 21 23 0 ./CppUtilities.cc:79 (parallel [
            (set (reg:SI 76)
                (and:SI (mem/s:SI (plus:DI (plus:DI (reg/f:DI 67 [ D.36005 ])
                                (reg:DI 75))
                            (const_int 32 [0x20])) [39 <variable>.D.25875._M_streambuf_state+0 S4 A64])
                    (const_int 5 [0x5])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 23 22 24 0 ./CppUtilities.cc:79 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 76)
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 24 23 27 0 ./CppUtilities.cc:79 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 74)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3300 [0xce4])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 27 24 25 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 25 27 26 1 ./CppUtilities.cc:79 (set (pc)
        (label_ref 59)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 26 25 28)

;; Start of basic block 2, registers live: (nil)
(code_label 28 26 29 2 99 "" [1 uses])

(note 29 28 30 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 30 29 31 2 ("./CppUtilities.cc") 85)

(insn 31 30 32 2 ./CppUtilities.cc:85 (parallel [
            (set (reg:DI 59 [ temp.612 ])
                (ashift:DI (reg:DI 60 [ temp.611 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 32 31 33 2 ./CppUtilities.cc:85 (set (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])
        (reg:DI 59 [ temp.612 ])) -1 (nil)
    (nil))

(note 33 32 34 2 ("./CppUtilities.cc") 86)

(insn 34 33 35 2 ./CppUtilities.cc:86 (set (reg:DI 4 si)
        (reg:DI 59 [ temp.612 ])) -1 (nil)
    (nil))

(insn 35 34 36 2 ./CppUtilities.cc:86 (set (reg:DI 5 di)
        (reg:DI 62 [ temp.607 ])) -1 (nil)
    (nil))

(call_insn 36 35 37 2 ./CppUtilities.cc:86 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Z11SafeReallocPvm") [flags 0x3] <function_decl 0x2b4c4365c300 SafeRealloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 37 36 38 2 ./CppUtilities.cc:86 (set (reg:DI 65 [ D.36023 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 38 37 39 2 ./CppUtilities.cc:86 (set (reg:DI 61 [ temp.608 ])
        (reg:DI 65 [ D.36023 ])) -1 (nil)
    (nil))

(insn 39 38 40 2 ./CppUtilities.cc:86 (set (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])
        (reg:DI 61 [ temp.608 ])) -1 (nil)
    (nil))

(note 40 39 41 2 ("./CppUtilities.cc") 88)

(insn 41 40 42 2 ./CppUtilities.cc:88 (set (reg:DI 77)
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 42 41 43 2 ./CppUtilities.cc:88 (set (mem:QI (plus:DI (plus:DI (reg:DI 61 [ temp.608 ])
                    (reg:DI 77))
                (const_int -1 [0xffffffffffffffff])) [0 S1 A8])
        (const_int 36 [0x24])) -1 (nil)
    (nil))

(note 43 42 44 2 ("./CppUtilities.cc") 89)

(insn 44 43 45 2 ./CppUtilities.cc:89 (set (reg:DI 64 [ old.99 ])
        (reg:DI 60 [ temp.611 ])) -1 (nil)
    (nil))

(insn 45 44 46 2 ./CppUtilities.cc:89 (parallel [
            (set (reg:DI 78)
                (plus:DI (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 46 45 47 2 ./CppUtilities.cc:89 (parallel [
            (set (reg:DI 79)
                (minus:DI (reg:DI 78)
                    (reg:DI 64 [ old.99 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 47 46 48 2 ./CppUtilities.cc:89 (parallel [
            (set (reg:DI 80)
                (plus:DI (reg:DI 64 [ old.99 ])
                    (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 48 47 49 2 ./CppUtilities.cc:89 (parallel [
            (set (reg:DI 81)
                (plus:DI (reg:DI 80)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 49 48 50 2 ./CppUtilities.cc:89 (set (reg:DI 1 dx)
        (reg:DI 79)) -1 (nil)
    (nil))

(insn 50 49 51 2 ./CppUtilities.cc:89 (set (reg:DI 4 si)
        (reg:DI 81)) -1 (nil)
    (nil))

(insn 51 50 52 2 ./CppUtilities.cc:89 (set (reg:DI 5 di)
        (reg/v/f:DI 69 [ file ])) -1 (nil)
    (nil))

(call_insn 52 51 53 2 ./CppUtilities.cc:89 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSi4readEPcl") [flags 0x41] <function_decl 0x2b4c433c7400 read>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 53 52 54 2 ./CppUtilities.cc:89 (set (reg/f:DI 63 [ D.36033 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 54 53 55 2 ./CppUtilities.cc:89 (set (reg/f:DI 82)
        (mem/s/f:DI (reg/f:DI 63 [ D.36033 ]) [6 <variable>._vptr.basic_istream+0 S8 A64])) -1 (nil)
    (nil))

(insn 55 54 56 2 ./CppUtilities.cc:89 (set (reg:DI 83)
        (mem:DI (plus:DI (reg/f:DI 82)
                (const_int -24 [0xffffffffffffffe8])) [19 S8 A64])) -1 (nil)
    (nil))

(insn 56 55 57 2 ./CppUtilities.cc:89 (parallel [
            (set (reg:SI 84)
                (and:SI (mem/s:SI (plus:DI (plus:DI (reg/f:DI 63 [ D.36033 ])
                                (reg:DI 83))
                            (const_int 32 [0x20])) [39 <variable>.D.25875._M_streambuf_state+0 S4 A64])
                    (const_int 5 [0x5])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 57 56 58 2 ./CppUtilities.cc:89 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 84)
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 58 57 59 2 ./CppUtilities.cc:89 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 74)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1100 [0x44c])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 59 58 60 3 98 "" [1 uses])

(note 60 59 61 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 61 60 62 3 ("./CppUtilities.cc") 82)

(insn 62 61 63 3 ./CppUtilities.cc:82 (set (reg:DI 62 [ temp.607 ])
        (mem/f:DI (reg/v/f:DI 70 [ buff ]) [11 S8 A64])) -1 (nil)
    (nil))

(insn 63 62 64 3 ./CppUtilities.cc:82 (set (reg:DI 60 [ temp.611 ])
        (mem:DI (reg/v/f:DI 71 [ size ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 64 63 65 3 ./CppUtilities.cc:82 (set (reg:DI 58 [ temp.617 ])
        (reg:DI 60 [ temp.611 ])) -1 (nil)
    (nil))

(insn 65 64 66 3 ./CppUtilities.cc:82 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:QI (plus:DI (plus:DI (reg:DI 62 [ temp.607 ])
                        (reg:DI 58 [ temp.617 ]))
                    (const_int -1 [0xffffffffffffffff])) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 66 65 68 3 ./CppUtilities.cc:82 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 79)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 68 66 69 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 69 68 70 4 ./CppUtilities.cc:82 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:QI (plus:DI (plus:DI (reg:DI 62 [ temp.607 ])
                        (reg:DI 58 [ temp.617 ]))
                    (const_int -2 [0xfffffffffffffffe])) [0 S1 A8])
            (const_int 10 [0xa]))) -1 (nil)
    (nil))

(jump_insn 70 69 73 4 ./CppUtilities.cc:82 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 28)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 73 70 71 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 71 73 72 5 ./CppUtilities.cc:82 (set (pc)
        (label_ref 79)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 72 71 74)

;; Start of basic block 6, registers live: (nil)
(code_label 74 72 75 6 97 "" [2 uses])

(note 75 74 76 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 76 75 77 6 ./CppUtilities.cc:82 (set (reg:SI 66 [ D.36014 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 77 76 78 6 ./CppUtilities.cc:82 (set (pc)
        (label_ref 82)) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

(barrier 78 77 79)

;; Start of basic block 7, registers live: (nil)
(code_label 79 78 80 7 100 "" [2 uses])

(note 80 79 81 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 81 80 82 7 ./CppUtilities.cc:82 (set (reg:SI 66 [ D.36014 ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 82 81 83 8 102 "" [1 uses])

(note 83 82 84 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 84 83 87 8 ./CppUtilities.cc:82 (set (reg:SI 68 [ <result> ])
        (reg:SI 66 [ D.36014 ])) -1 (nil)
    (nil))

(note 87 84 88 8 NOTE_INSN_FUNCTION_END)

(note 88 87 90 8 ("./CppUtilities.cc") 94)

(insn 90 88 96 8 ./CppUtilities.cc:94 (set (reg/i:SI 0 ax)
        (reg:SI 68 [ <result> ])) -1 (nil)
    (nil))

(insn 96 90 0 8 ./CppUtilities.cc:94 (use (reg/i:SI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)


;; Function std::string pathWithoutSuffix(const std::string&) (_Z17pathWithoutSuffixRKSs)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Merged 2 and 3 without moving.
Edge 2->31 redirected to 32
Deleted label in block 4.
Deleted label in block 5.
Merged 6 and 7 without moving.
Forwarding edge 6->8 to 32 failed.
Redirecting jump 49 from 31 to 32.
Deleted label in block 10.
Merged 11 and 12 without moving.
Deleted label in block 13.
Merged 18 and 19 without moving.
Edge 18->31 redirected to 32
Deleted label in block 20.
Redirecting jump 123 from 31 to 32.
Deleted label in block 24.
Merged 25 and 26 without moving.
Deleted label in block 27.
Deleting block 31.
Merged 32 and 33 without moving.


try_optimize_cfg iteration 2

Forwarding edge 6->8 to 32 failed.


try_optimize_cfg iteration 1

Forwarding edge 4->5 to 30 failed.
Merged 6 and 7 without moving.
Merged 11 and 12 without moving.
Merged 18 and 19 without moving.
Merged 18 and 20 without moving.
Merged 26 and 27 without moving.


try_optimize_cfg iteration 2

Forwarding edge 4->5 to 30 failed.
(note 1 0 7 ("./CppUtilities.cc") 184)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppUtilities.cc:184 (set (reg/f:DI 86 [ <result> ])
        (reg:DI 5 di [ D.39584 ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:184 (set (reg/v/f:DI 87 [ file ])
        (reg:DI 4 si [ file ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppUtilities.cc") 186)

(insn 10 9 11 0 ./CppUtilities.cc:186 (set (reg/f:DI 85 [ base.103 ])
        (reg/f:DI 86 [ <result> ])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppUtilities.cc:186 (set (reg:DI 4 si)
        (reg/v/f:DI 87 [ file ])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppUtilities.cc:186 (set (reg:DI 5 di)
        (reg/f:DI 85 [ base.103 ])) -1 (nil)
    (nil))

(call_insn 13 12 14 0 ./CppUtilities.cc:186 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSs") [flags 0x41] <function_decl 0x2b4c42b2a900 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 14 13 15 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1773)

(insn 15 14 16 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 1 dx)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 16 15 17 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:SI 4 si)
        (const_int 46 [0x2e])) -1 (nil)
    (nil))

(insn 17 16 18 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 5 di)
        (reg/v/f:DI 87 [ file ])) -1 (nil)
    (nil))

(call_insn 18 17 210 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNKSs5rfindEcm") [flags 0x41] <function_decl 0x2b4c42b13500 rfind>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 210 18 19 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 19 210 21 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg/v:DI 84 [ ind ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 21 19 22 1 ("./CppUtilities.cc") 189)

(insn 22 21 23 1 ./CppUtilities.cc:189 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 84 [ ind ])
            (const_int -1 [0xffffffffffffffff]))) -1 (nil)
    (nil))

(jump_insn 23 22 25 1 ./CppUtilities.cc:189 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 200)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 3195 [0xc7b])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 25 23 26 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 26 25 27 2 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1908)

(insn 27 26 28 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 88)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -16 [0xfffffffffffffff0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 28 27 29 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 2 cx)
        (reg/v:DI 84 [ ind ])) -1 (nil)
    (nil))

(insn 29 28 30 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 1 dx)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 30 29 31 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 4 si)
        (reg/v/f:DI 87 [ file ])) -1 (nil)
    (nil))

(insn 31 30 32 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 5 di)
        (reg:DI 88)) -1 (nil)
    (nil))

(call_insn 32 31 34 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSsmm") [flags 0x41] <function_decl 0x2b4c42b2a700 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 34 32 35 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 35 34 36 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 36 35 37 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 89)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -16 [0xfffffffffffffff0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 37 36 38 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 89)) -1 (nil)
    (nil))

(insn 38 37 39 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 85 [ base.103 ])) -1 (nil)
    (nil))

(call_insn 39 38 211 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 2 [0x2])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 211 39 40 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 40 211 42 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 79 [ D.39366 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 42 40 43 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 43 42 44 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 65 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -16 [0xfffffffffffffff0])) [11 D.36176._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 44 43 45 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 45 44 46 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 66 [ D.39468 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 46 45 47 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 47 46 48 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 65 [ this ])
            (reg:DI 66 [ D.39468 ]))) -1 (nil)
    (nil))

(jump_insn 48 47 51 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 105)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 99 [0x63])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 99 [0x63]))
            (nil))))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 51 48 49 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 49 51 50 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (label_ref:DI 200)) 548 {jump} (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 50 49 259)

;; Start of basic block 6, registers live: (nil)
(code_label/s 259 50 262 6 133 "" [1 uses])

(note 262 259 260 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 260 262 261 6 (set (reg:DI 91)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 261 260 52 6 (set (reg:DI 90)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 52 261 54 6 "" NOTE_INSN_DELETED_LABEL 111)

(insn 54 52 55 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:SI 82 [ save_filt.156 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 55 54 56 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:DI 83 [ save_eptr.155 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(note 56 55 57 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 57 56 58 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 72 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -16 [0xfffffffffffffff0])) [11 D.36176._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 58 57 59 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 59 58 60 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 73 [ D.39412 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 60 59 61 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 61 60 62 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 72 [ this ])
            (reg:DI 73 [ D.39412 ]))) -1 (nil)
    (nil))

(jump_insn 62 61 64 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 91)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 64 62 65 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 65 64 66 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 66 65 67 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 92)
                (plus:DI (reg:DI 72 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 67 66 68 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 68 67 69 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 92)) -1 (nil)
    (nil))

(call_insn 69 68 212 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 212 69 70 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 70 212 72 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 74 [ D.39420 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 72 70 73 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 74 [ D.39420 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 73 72 75 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 91)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(note 75 73 76 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 76 75 77 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 77 76 78 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 93)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -3 [0xfffffffffffffffd])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 78 77 79 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 93)) -1 (nil)
    (nil))

(insn 79 78 80 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 72 [ this ])) -1 (nil)
    (nil))

(call_insn 80 79 81 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 81 80 82 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 91)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 82 81 251)

;; Start of basic block 11, registers live: (nil)
(code_label/s 251 82 254 11 131 "" [1 uses])

(note 254 251 252 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 252 254 253 11 (set (reg:DI 91)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 253 252 83 11 (set (reg:DI 90)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 83 253 85 11 "" NOTE_INSN_DELETED_LABEL 115)

(insn 85 83 86 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 76 [ save_filt.160 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 86 85 87 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 75 [ save_eptr.159 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(insn 87 86 88 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 75 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 88 87 226 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 76 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 226 88 90 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 97)) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

(barrier 90 226 91)

;; Start of basic block 13, registers live: (nil)
(code_label 91 90 92 13 112 "" [3 uses])

(note 92 91 93 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 93 92 94 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 83 [ save_eptr.155 ])) -1 (nil)
    (nil))

(insn 94 93 230 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 82 [ save_filt.156 ]))) -1 (nil)
    (nil))

(jump_insn 230 94 96 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 141)) -1 (nil)
    (nil))
;; End of basic block 13, registers live:
 (nil)

(barrier 96 230 255)

(note/s 255 96 97 39 "" NOTE_INSN_DELETED_LABEL 132)

;; Start of basic block 14, registers live: (nil)
(code_label/s 97 255 98 14 116 "" [2 uses])

(note 98 97 99 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(insn 99 98 100 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 78 [ save_filt.162 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 100 99 101 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 77 [ save_eptr.161 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(insn 101 100 102 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 77 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 102 101 228 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 78 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 228 102 104 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 141)) -1 (nil)
    (nil))
;; End of basic block 14, registers live:
 (nil)

(barrier 104 228 105)

;; Start of basic block 15, registers live: (nil)
(code_label 105 104 106 15 110 "" [1 uses])

(note 106 105 107 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 107 106 108 15 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 108 107 109 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 94)
                (plus:DI (reg:DI 65 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 109 108 110 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 110 109 111 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 94)) -1 (nil)
    (nil))

(call_insn 111 110 213 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 10 [0xa])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(note 213 111 112 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(insn 112 213 114 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 67 [ D.39476 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 114 112 115 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 67 [ D.39476 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 115 114 117 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 200)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 8122 [0x1fba])
        (nil)))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(note 117 115 118 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(note 118 117 119 17 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 119 118 120 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 95)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -2 [0xfffffffffffffffe])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 120 119 121 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 95)) -1 (nil)
    (nil))

(insn 121 120 122 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 65 [ this ])) -1 (nil)
    (nil))

(call_insn 122 121 123 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 123 122 124 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref:DI 200)) 548 {jump} (nil)
    (nil))
;; End of basic block 17, registers live:
 (nil)

(barrier 124 123 243)

;; Start of basic block 18, registers live: (nil)
(code_label/s 243 124 246 18 129 "" [1 uses])

(note 246 243 244 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(insn 244 246 245 18 (set (reg:DI 91)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 245 244 125 18 (set (reg:DI 90)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 125 245 127 18 "" NOTE_INSN_DELETED_LABEL 118)

(insn 127 125 128 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 69 [ save_filt.160 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 128 127 129 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 68 [ save_eptr.159 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(insn 129 128 130 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 68 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 130 129 247 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 69 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 247 130 133 18 "" NOTE_INSN_DELETED_LABEL 130)

(note/s 133 247 135 18 "" NOTE_INSN_DELETED_LABEL 119)

(insn 135 133 136 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 71 [ save_filt.162 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 136 135 137 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 70 [ save_eptr.161 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(insn 137 136 138 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 70 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 138 137 224 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 71 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 224 138 140 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 141)) -1 (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

(barrier 140 224 263)

;; Start of basic block 21, registers live: (nil)
(code_label/s 263 140 266 21 134 "" [1 uses])

(note 266 263 264 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(insn 264 266 265 21 (set (reg:DI 91)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 265 264 141 21 (set (reg:DI 90)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 21, registers live:
 (nil)

;; Start of basic block 22, registers live: (nil)
(code_label/s 141 265 142 22 120 "" [4 uses])

(note 142 141 143 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(insn 143 142 144 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 80 [ save_filt.158 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 144 143 145 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 81 [ save_eptr.157 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(note 145 144 146 22 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 146 145 147 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 58 [ this ])
                (plus:DI (mem/s/f:DI (reg/f:DI 85 [ base.103 ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 147 146 148 22 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 148 147 149 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 59 [ D.39525 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 149 148 150 22 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 150 149 151 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 58 [ this ])
            (reg:DI 59 [ D.39525 ]))) -1 (nil)
    (nil))

(jump_insn 151 150 153 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 180)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 22, registers live:
 (nil)

;; Start of basic block 23, registers live: (nil)
(note 153 151 154 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(note 154 153 155 23 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 155 154 156 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 96)
                (plus:DI (reg:DI 58 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 156 155 157 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 157 156 158 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 96)) -1 (nil)
    (nil))

(call_insn 158 157 214 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 14 [0xe])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 23, registers live:
 (nil)

;; Start of basic block 24, registers live: (nil)
(note 214 158 159 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(insn 159 214 161 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 60 [ D.39533 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 161 159 162 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 60 [ D.39533 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 162 161 164 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 180)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(note 164 162 165 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(note 165 164 166 25 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 166 165 167 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 97)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 167 166 168 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 97)) -1 (nil)
    (nil))

(insn 168 167 169 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 58 [ this ])) -1 (nil)
    (nil))

(call_insn 169 168 170 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 170 169 171 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 180)) -1 (nil)
    (nil))
;; End of basic block 25, registers live:
 (nil)

(barrier 171 170 235)

;; Start of basic block 26, registers live: (nil)
(code_label/s 235 171 238 26 127 "" [1 uses])

(note 238 235 236 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(insn 236 238 237 26 (set (reg:DI 91)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 237 236 172 26 (set (reg:DI 90)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 172 237 174 26 "" NOTE_INSN_DELETED_LABEL 124)

(insn 174 172 175 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 62 [ save_filt.160 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 175 174 176 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 61 [ save_eptr.159 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(insn 176 175 177 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 61 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 177 176 217 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 62 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 217 177 179 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 186)) -1 (nil)
    (nil))
;; End of basic block 26, registers live:
 (nil)

(barrier 179 217 180)

;; Start of basic block 28, registers live: (nil)
(code_label 180 179 181 28 121 "" [3 uses])

(note 181 180 182 28 [bb 28] NOTE_INSN_BASIC_BLOCK)

(insn 182 181 183 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 81 [ save_eptr.157 ])) -1 (nil)
    (nil))

(insn 183 182 232 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 80 [ save_filt.158 ]))) -1 (nil)
    (nil))

(insn 232 183 233 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 91)) -1 (nil)
    (nil))

(call_insn 233 232 185 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 28, registers live:
 (nil)

(barrier 185 233 239)

(note/s 239 185 186 35 "" NOTE_INSN_DELETED_LABEL 128)

;; Start of basic block 29, registers live: (nil)
(code_label/s 186 239 187 29 125 "" [2 uses])

(note 187 186 188 29 [bb 29] NOTE_INSN_BASIC_BLOCK)

(insn 188 187 189 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 64 [ save_filt.162 ])
        (subreg:SI (reg:DI 90) 0)) -1 (nil)
    (nil))

(insn 189 188 190 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 63 [ save_eptr.161 ])
        (reg:DI 91)) -1 (nil)
    (nil))

(insn 190 189 191 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 91)
        (reg:DI 63 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 191 190 219 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90)
        (sign_extend:DI (reg:SI 64 [ save_filt.162 ]))) -1 (nil)
    (nil))

(insn 219 191 220 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 91)) -1 (nil)
    (nil))

(call_insn 220 219 193 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 29, registers live:
 (nil)

(barrier 193 220 198)

(note 198 193 199 NOTE_INSN_FUNCTION_END)

(note 199 198 200 ("./CppUtilities.cc") 193)

;; Start of basic block 30, registers live: (nil)
(code_label 200 199 215 30 106 "" [4 uses])

(note 215 200 201 30 [bb 30] NOTE_INSN_BASIC_BLOCK)

(insn 201 215 202 30 ./CppUtilities.cc:193 (set (reg/i:DI 0 ax)
        (reg/f:DI 86 [ <result> ])) -1 (nil)
    (nil))

(insn 202 201 208 30 ./CppUtilities.cc:193 (set (reg/i:DI 0 ax)
        (reg/f:DI 86 [ <result> ])) -1 (nil)
    (nil))

(insn 208 202 0 30 ./CppUtilities.cc:193 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 30, registers live:
 (nil)


;; Function std::string baseFileName(const std::string&, bool) (_Z12baseFileNameRKSsb)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Merged 2 and 3 without moving.
Deleted label in block 4.
Deleted label in block 5.
Deleted label in block 6.
Deleted label in block 8.
Merged 9 and 10 without moving.
Edge 9->99 redirected to 100
Deleted label in block 11.
Merged 12 and 13 without moving.
Edge 12->99 redirected to 100
Deleted label in block 14.
Redirecting jump 87 from 99 to 100.
Deleted label in block 18.
Merged 19 and 20 without moving.
Deleted label in block 21.
Merged 26 and 27 without moving.
Deleted label in block 29.
Deleted label in block 30.
Deleted label in block 32.
Merged 33 and 34 without moving.
Edge 33->99 redirected to 100
Deleted label in block 35.
Merged 36 and 37 without moving.
Edge 36->99 redirected to 100
Deleted label in block 38.
Redirecting jump 237 from 99 to 100.
Deleted label in block 42.
Merged 43 and 44 without moving.
Deleted label in block 45.
Deleted label in block 51.
Deleted label in block 52.
Merged 53 and 54 without moving.
Edge 53->99 redirected to 100
Deleted label in block 55.
Merged 56 and 57 without moving.
Edge 56->99 redirected to 100
Deleted label in block 58.
Redirecting jump 358 from 99 to 100.
Deleted label in block 62.
Merged 63 and 64 without moving.
Deleted label in block 65.
Edge 69->99 redirected to 100
Edge 70->99 redirected to 100
Deleted label in block 71.
Deleted label in block 72.
Deleted label in block 74.
Merged 75 and 76 without moving.
Edge 75->99 redirected to 100
Deleted label in block 77.
Merged 78 and 79 without moving.
Edge 78->99 redirected to 100
Deleted label in block 80.
Redirecting jump 493 from 99 to 100.
Deleted label in block 84.
Merged 85 and 86 without moving.
Deleted label in block 87.
Deleted label in block 92.
Merged 93 and 94 without moving.
Deleted label in block 95.
Deleting block 99.
Merged 100 and 101 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

Merged 11 and 12 without moving.
Merged 11 and 13 without moving.
Merged 14 and 15 without moving.
Merged 19 and 20 without moving.
Merged 34 and 35 without moving.
Merged 34 and 36 without moving.
Merged 37 and 38 without moving.
Merged 42 and 43 without moving.
Merged 54 and 55 without moving.
Merged 54 and 56 without moving.
Merged 57 and 58 without moving.
Merged 62 and 63 without moving.
Merged 76 and 77 without moving.
Merged 76 and 78 without moving.
Merged 79 and 80 without moving.
Merged 84 and 85 without moving.
Merged 93 and 94 without moving.


try_optimize_cfg iteration 2

(note 1 0 9 ("./CppUtilities.cc") 148)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 ./CppUtilities.cc:148 (set (reg/f:DI 143 [ <result> ])
        (reg:DI 5 di [ D.40344 ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:148 (set (reg/v/f:DI 144 [ file ])
        (reg:DI 4 si [ file ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppUtilities.cc:148 (set (reg:SI 146)
        (reg:SI 1 dx [ withSuffix ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppUtilities.cc:148 (set (reg/v:QI 145 [ withSuffix ])
        (subreg:QI (reg:SI 146) 0)) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("./CppUtilities.cc") 150)

(insn 12 11 13 0 ./CppUtilities.cc:150 (set (reg/f:DI 142 [ base.102 ])
        (reg/f:DI 143 [ <result> ])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppUtilities.cc:150 (set (reg:DI 4 si)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppUtilities.cc:150 (set (reg:DI 5 di)
        (reg/f:DI 142 [ base.102 ])) -1 (nil)
    (nil))

(call_insn 15 14 16 0 ./CppUtilities.cc:150 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSs") [flags 0x41] <function_decl 0x2b4c42b2a900 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(note 16 15 17 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1773)

(insn 17 16 18 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 1 dx)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 18 17 19 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:SI 4 si)
        (const_int 47 [0x2f])) -1 (nil)
    (nil))

(insn 19 18 20 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 5 di)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(call_insn 20 19 633 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNKSs5rfindEcm") [flags 0x41] <function_decl 0x2b4c42b13500 rfind>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 633 20 21 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 21 633 23 1 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg/v:DI 139 [ ind1 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 23 21 24 1 ("./CppUtilities.cc") 154)

(insn 24 23 25 1 ./CppUtilities.cc:154 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:QI 145 [ withSuffix ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 25 24 27 1 ./CppUtilities.cc:154 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 158)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 3000 [0xbb8])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 27 25 28 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 28 27 29 2 ./CppUtilities.cc:154 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 139 [ ind1 ])
            (const_int -1 [0xffffffffffffffff]))) -1 (nil)
    (nil))

(jump_insn 29 28 31 2 ./CppUtilities.cc:154 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 158)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 2900 [0xb54])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 31 29 32 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 32 31 33 3 ("./CppUtilities.cc") 156)

(insn 33 32 34 3 ./CppUtilities.cc:156 (parallel [
            (set (reg/v:DI 123 [ __pos ])
                (plus:DI (reg/v:DI 139 [ ind1 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 34 33 35 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 305)

(insn 35 34 36 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (reg/f:DI 147)
        (mem/s/f:DI (reg/v/f:DI 144 [ file ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil)
    (nil))

(insn 36 35 37 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 123 [ __pos ])
            (mem/s:DI (plus:DI (reg/f:DI 147)
                    (const_int -24 [0xffffffffffffffe8])) [19 <variable>.D.17834._M_length+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 37 36 39 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 44)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 39 37 40 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 40 39 41 4 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 306)

(insn 41 40 42 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:306 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC7") [flags 0x2] <string_cst 0x2b4c43cf0380>)) -1 (nil)
    (nil))

(call_insn 42 41 43 4 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:306 (call (mem:QI (symbol_ref:DI ("_ZSt20__throw_out_of_rangePKc") [flags 0x41] <function_decl 0x2b4c423ce000 __throw_out_of_range>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 1 [0x1])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 4, registers live:
 (nil)

(barrier 43 42 44)

;; Start of basic block 5, registers live: (nil)
(code_label 44 43 45 5 142 "" [1 uses])

(note 45 44 46 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 46 45 47 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1908)

(insn 47 46 48 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 148)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -80 [0xffffffffffffffb0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 48 47 49 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 2 cx)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 49 48 50 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 1 dx)
        (reg/v:DI 123 [ __pos ])) -1 (nil)
    (nil))

(insn 50 49 51 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 4 si)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(insn 51 50 52 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 5 di)
        (reg:DI 148)) -1 (nil)
    (nil))

(call_insn 52 51 54 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSsmm") [flags 0x41] <function_decl 0x2b4c42b2a700 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(note 54 52 55 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 55 54 56 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 56 55 57 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 149)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -80 [0xffffffffffffffb0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 57 56 58 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 149)) -1 (nil)
    (nil))

(insn 58 57 59 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 142 [ base.102 ])) -1 (nil)
    (nil))

(call_insn 59 58 634 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 2 [0x2])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 634 59 60 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 60 634 62 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 127 [ D.38138 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 62 60 63 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 63 62 64 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 116 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -80 [0xffffffffffffffb0])) [11 D.36147._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 64 63 65 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 65 64 66 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 117 [ D.39737 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 66 65 67 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 67 66 68 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 116 [ this ])
            (reg:DI 117 [ D.39737 ]))) -1 (nil)
    (nil))

(jump_insn 68 67 70 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 70 68 71 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 71 70 72 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 72 71 73 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 150)
                (plus:DI (reg:DI 116 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 73 72 74 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 74 73 75 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 150)) -1 (nil)
    (nil))

(call_insn 75 74 635 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 12 [0xc])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 635 75 76 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 76 635 78 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 118 [ D.39745 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 78 76 79 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 118 [ D.39745 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 79 78 81 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 8122 [0x1fba])
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(note 81 79 82 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 82 81 83 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 83 82 84 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 151)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -9 [0xfffffffffffffff7])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 84 83 85 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 151)) -1 (nil)
    (nil))

(insn 85 84 86 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 116 [ this ])) -1 (nil)
    (nil))

(call_insn 86 85 87 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 87 86 88 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref:DI 623)) 548 {jump} (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 88 87 765)

;; Start of basic block 11, registers live: (nil)
(code_label/s 765 88 768 11 221 "" [1 uses])

(note 768 765 766 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 766 768 767 11 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 767 766 89 11 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 89 767 91 11 "" NOTE_INSN_DELETED_LABEL 148)

(insn 91 89 92 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 120 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 92 91 93 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 119 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 93 92 94 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 119 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 94 93 769 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 120 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 769 94 97 11 "" NOTE_INSN_DELETED_LABEL 222)

(note/s 97 769 99 11 "" NOTE_INSN_DELETED_LABEL 149)

(insn 99 97 100 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 122 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 100 99 101 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 121 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 101 100 102 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 121 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 102 101 688 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 122 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 688 102 104 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

(barrier 104 688 785)

;; Start of basic block 14, registers live: (nil)
(code_label/s 785 104 788 14 226 "" [1 uses])

(note 788 785 786 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(insn 786 788 787 14 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 787 786 105 14 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 105 787 107 14 "" NOTE_INSN_DELETED_LABEL 150)

(insn 107 105 108 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 136 [ save_filt.164 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 108 107 109 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 137 [ save_eptr.163 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(note 109 108 110 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 110 109 111 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 109 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -80 [0xffffffffffffffb0])) [11 D.36147._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 111 110 112 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 112 111 113 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 110 [ D.39793 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 113 112 114 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 114 113 115 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 109 [ this ])
            (reg:DI 110 [ D.39793 ]))) -1 (nil)
    (nil))

(jump_insn 115 114 117 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 144)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 14, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(note 117 115 118 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(note 118 117 119 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 119 118 120 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 154)
                (plus:DI (reg:DI 109 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 120 119 121 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 121 120 122 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 154)) -1 (nil)
    (nil))

(call_insn 122 121 636 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 16 [0x10])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(note 636 122 123 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(insn 123 636 125 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 111 [ D.39801 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 125 123 126 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 111 [ D.39801 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 126 125 128 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 144)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(note 128 126 129 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(note 129 128 130 18 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 130 129 131 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 155)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -8 [0xfffffffffffffff8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 131 130 132 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 155)) -1 (nil)
    (nil))

(insn 132 131 133 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 109 [ this ])) -1 (nil)
    (nil))

(call_insn 133 132 134 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 134 133 135 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 144)) -1 (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

(barrier 135 134 757)

;; Start of basic block 19, registers live: (nil)
(code_label/s 757 135 760 19 219 "" [1 uses])

(note 760 757 758 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(insn 758 760 759 19 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 759 758 136 19 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 136 759 138 19 "" NOTE_INSN_DELETED_LABEL 154)

(insn 138 136 139 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 113 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 139 138 140 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 112 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 140 139 141 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 112 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 141 140 682 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 113 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 682 141 143 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 150)) -1 (nil)
    (nil))
;; End of basic block 19, registers live:
 (nil)

(barrier 143 682 144)

;; Start of basic block 21, registers live: (nil)
(code_label 144 143 145 21 151 "" [3 uses])

(note 145 144 146 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(insn 146 145 147 21 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 137 [ save_eptr.163 ])) -1 (nil)
    (nil))

(insn 147 146 696 21 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 136 [ save_filt.164 ]))) -1 (nil)
    (nil))

(jump_insn 696 147 149 21 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 21, registers live:
 (nil)

(barrier 149 696 761)

(note/s 761 149 150 117 "" NOTE_INSN_DELETED_LABEL 220)

;; Start of basic block 22, registers live: (nil)
(code_label/s 150 761 151 22 155 "" [2 uses])

(note 151 150 152 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(insn 152 151 153 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 115 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 153 152 154 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 114 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 154 153 155 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 114 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 155 154 684 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 115 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 684 155 157 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 22, registers live:
 (nil)

(barrier 157 684 158)

;; Start of basic block 23, registers live: (nil)
(code_label 158 157 159 23 139 "" [2 uses])

(note 159 158 160 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(note 160 159 161 23 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1773)

(insn 161 160 162 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 1 dx)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 162 161 163 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:SI 4 si)
        (const_int 46 [0x2e])) -1 (nil)
    (nil))

(insn 163 162 164 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 5 di)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(call_insn 164 163 637 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNKSs5rfindEcm") [flags 0x41] <function_decl 0x2b4c42b13500 rfind>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 23, registers live:
 (nil)

;; Start of basic block 24, registers live: (nil)
(note 637 164 165 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(insn 165 637 167 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1773 (set (reg/v:DI 138 [ ind2 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 167 165 168 24 ("./CppUtilities.cc") 163)

(insn 168 167 169 24 ./CppUtilities.cc:163 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 139 [ ind1 ])
            (const_int -1 [0xffffffffffffffff]))) -1 (nil)
    (nil))

(insn 169 168 170 24 ./CppUtilities.cc:163 (set (reg:QI 156)
        (ne:QI (reg:CCZ 17 flags)
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EQUAL (ne:QI (reg/v:DI 139 [ ind1 ])
            (const_int -1 [0xffffffffffffffff]))
        (nil)))

(insn 170 169 171 24 ./CppUtilities.cc:163 (set (reg:QI 141 [ D.36161 ])
        (reg:QI 156)) -1 (nil)
    (nil))

(insn 171 170 172 24 ./CppUtilities.cc:163 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 138 [ ind2 ])
            (const_int -1 [0xffffffffffffffff]))) -1 (nil)
    (nil))

(insn 172 171 173 24 ./CppUtilities.cc:163 (set (reg:QI 157)
        (ne:QI (reg:CCZ 17 flags)
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EQUAL (ne:QI (reg/v:DI 138 [ ind2 ])
            (const_int -1 [0xffffffffffffffff]))
        (nil)))

(insn 173 172 174 24 ./CppUtilities.cc:163 (set (reg:QI 140 [ D.36162 ])
        (reg:QI 157)) -1 (nil)
    (nil))

(insn 174 173 175 24 ./CppUtilities.cc:163 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 141 [ D.36161 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 175 174 638 24 ./CppUtilities.cc:163 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 308)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(note 638 175 176 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(insn 176 638 177 25 ./CppUtilities.cc:163 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 140 [ D.36162 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 177 176 179 25 ./CppUtilities.cc:163 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 308)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 25, registers live:
 (nil)

;; Start of basic block 26, registers live: (nil)
(note 179 177 180 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(note 180 179 181 26 ("./CppUtilities.cc") 165)

(insn 181 180 182 26 ./CppUtilities.cc:165 (parallel [
            (set (reg/v:DI 108 [ __pos ])
                (plus:DI (reg/v:DI 139 [ ind1 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 182 181 183 26 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 305)

(insn 183 182 184 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (reg/f:DI 158)
        (mem/s/f:DI (reg/v/f:DI 144 [ file ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil)
    (nil))

(insn 184 183 185 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 108 [ __pos ])
            (mem/s:DI (plus:DI (reg/f:DI 158)
                    (const_int -24 [0xffffffffffffffe8])) [19 <variable>.D.17834._M_length+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 185 184 187 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 192)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 26, registers live:
 (nil)

;; Start of basic block 27, registers live: (nil)
(note 187 185 188 27 [bb 27] NOTE_INSN_BASIC_BLOCK)

(note 188 187 189 27 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 306)

(insn 189 188 190 27 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:306 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC7") [flags 0x2] <string_cst 0x2b4c43cf0380>)) -1 (nil)
    (nil))

(call_insn 190 189 191 27 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:306 (call (mem:QI (symbol_ref:DI ("_ZSt20__throw_out_of_rangePKc") [flags 0x41] <function_decl 0x2b4c423ce000 __throw_out_of_range>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 1 [0x1])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 27, registers live:
 (nil)

(barrier 191 190 192)

;; Start of basic block 28, registers live: (nil)
(code_label 192 191 193 28 158 "" [1 uses])

(note 193 192 194 28 [bb 28] NOTE_INSN_BASIC_BLOCK)

(note 194 193 195 28 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1908)

(insn 195 194 196 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 159)
                (minus:DI (reg/v:DI 138 [ ind2 ])
                    (reg/v:DI 139 [ ind1 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 196 195 197 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 160)
                (plus:DI (reg:DI 159)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 197 196 198 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 161)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 198 197 199 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 2 cx)
        (reg:DI 160)) -1 (nil)
    (nil))

(insn 199 198 200 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 1 dx)
        (reg/v:DI 108 [ __pos ])) -1 (nil)
    (nil))

(insn 200 199 201 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 4 si)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(insn 201 200 202 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 5 di)
        (reg:DI 161)) -1 (nil)
    (nil))

(call_insn 202 201 204 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSsmm") [flags 0x41] <function_decl 0x2b4c42b2a700 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 28, registers live:
 (nil)

;; Start of basic block 29, registers live: (nil)
(note 204 202 205 29 [bb 29] NOTE_INSN_BASIC_BLOCK)

(note 205 204 206 29 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 206 205 207 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 162)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -64 [0xffffffffffffffc0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 207 206 208 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 162)) -1 (nil)
    (nil))

(insn 208 207 209 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 142 [ base.102 ])) -1 (nil)
    (nil))

(call_insn 209 208 639 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 4 [0x4])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 29, registers live:
 (nil)

;; Start of basic block 30, registers live: (nil)
(note 639 209 210 30 [bb 30] NOTE_INSN_BASIC_BLOCK)

(insn 210 639 212 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 126 [ D.38149 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 212 210 213 30 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 213 212 214 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 101 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -64 [0xffffffffffffffc0])) [11 D.36149._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 214 213 215 30 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 215 214 216 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 102 [ D.39871 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 216 215 217 30 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 217 216 218 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 101 [ this ])
            (reg:DI 102 [ D.39871 ]))) -1 (nil)
    (nil))

(jump_insn 218 217 220 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 30, registers live:
 (nil)

;; Start of basic block 31, registers live: (nil)
(note 220 218 221 31 [bb 31] NOTE_INSN_BASIC_BLOCK)

(note 221 220 222 31 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 222 221 223 31 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 163)
                (plus:DI (reg:DI 101 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 223 222 224 31 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 224 223 225 31 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 163)) -1 (nil)
    (nil))

(call_insn 225 224 640 31 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 20 [0x14])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 31, registers live:
 (nil)

;; Start of basic block 32, registers live: (nil)
(note 640 225 226 32 [bb 32] NOTE_INSN_BASIC_BLOCK)

(insn 226 640 228 32 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 103 [ D.39879 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 228 226 229 32 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 103 [ D.39879 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 229 228 231 32 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 8122 [0x1fba])
        (nil)))
;; End of basic block 32, registers live:
 (nil)

;; Start of basic block 33, registers live: (nil)
(note 231 229 232 33 [bb 33] NOTE_INSN_BASIC_BLOCK)

(note 232 231 233 33 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 233 232 234 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 164)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -7 [0xfffffffffffffff9])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 234 233 235 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 164)) -1 (nil)
    (nil))

(insn 235 234 236 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 101 [ this ])) -1 (nil)
    (nil))

(call_insn 236 235 237 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 237 236 238 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref:DI 623)) 548 {jump} (nil)
    (nil))
;; End of basic block 33, registers live:
 (nil)

(barrier 238 237 749)

;; Start of basic block 34, registers live: (nil)
(code_label/s 749 238 752 34 217 "" [1 uses])

(note 752 749 750 34 [bb 34] NOTE_INSN_BASIC_BLOCK)

(insn 750 752 751 34 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 751 750 239 34 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 239 751 241 34 "" NOTE_INSN_DELETED_LABEL 163)

(insn 241 239 242 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 105 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 242 241 243 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 104 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 243 242 244 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 104 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 244 243 753 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 105 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 753 244 247 34 "" NOTE_INSN_DELETED_LABEL 218)

(note/s 247 753 249 34 "" NOTE_INSN_DELETED_LABEL 164)

(insn 249 247 250 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 107 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 250 249 251 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 106 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 251 250 252 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 106 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 252 251 680 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 107 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 680 252 254 34 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 34, registers live:
 (nil)

(barrier 254 680 781)

;; Start of basic block 37, registers live: (nil)
(code_label/s 781 254 784 37 225 "" [1 uses])

(note 784 781 782 37 [bb 37] NOTE_INSN_BASIC_BLOCK)

(insn 782 784 783 37 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 783 782 255 37 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 255 783 257 37 "" NOTE_INSN_DELETED_LABEL 165)

(insn 257 255 258 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 134 [ save_filt.166 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 258 257 259 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 135 [ save_eptr.165 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(note 259 258 260 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 260 259 261 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 94 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -64 [0xffffffffffffffc0])) [11 D.36149._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 261 260 262 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 262 261 263 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 95 [ D.39927 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 263 262 264 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 264 263 265 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 94 [ this ])
            (reg:DI 95 [ D.39927 ]))) -1 (nil)
    (nil))

(jump_insn 265 264 267 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 294)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 37, registers live:
 (nil)

;; Start of basic block 39, registers live: (nil)
(note 267 265 268 39 [bb 39] NOTE_INSN_BASIC_BLOCK)

(note 268 267 269 39 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 269 268 270 39 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 165)
                (plus:DI (reg:DI 94 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 270 269 271 39 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 271 270 272 39 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 165)) -1 (nil)
    (nil))

(call_insn 272 271 641 39 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 24 [0x18])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 39, registers live:
 (nil)

;; Start of basic block 40, registers live: (nil)
(note 641 272 273 40 [bb 40] NOTE_INSN_BASIC_BLOCK)

(insn 273 641 275 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 96 [ D.39935 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 275 273 276 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 96 [ D.39935 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 276 275 278 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 294)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 40, registers live:
 (nil)

;; Start of basic block 41, registers live: (nil)
(note 278 276 279 41 [bb 41] NOTE_INSN_BASIC_BLOCK)

(note 279 278 280 41 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 280 279 281 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 166)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -6 [0xfffffffffffffffa])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 281 280 282 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 166)) -1 (nil)
    (nil))

(insn 282 281 283 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 94 [ this ])) -1 (nil)
    (nil))

(call_insn 283 282 284 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 284 283 285 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 294)) -1 (nil)
    (nil))
;; End of basic block 41, registers live:
 (nil)

(barrier 285 284 741)

;; Start of basic block 42, registers live: (nil)
(code_label/s 741 285 744 42 215 "" [1 uses])

(note 744 741 742 42 [bb 42] NOTE_INSN_BASIC_BLOCK)

(insn 742 744 743 42 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 743 742 286 42 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 286 743 288 42 "" NOTE_INSN_DELETED_LABEL 169)

(insn 288 286 289 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 98 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 289 288 290 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 97 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 290 289 291 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 97 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 291 290 674 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 98 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 674 291 293 42 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 300)) -1 (nil)
    (nil))
;; End of basic block 42, registers live:
 (nil)

(barrier 293 674 294)

;; Start of basic block 44, registers live: (nil)
(code_label 294 293 295 44 166 "" [3 uses])

(note 295 294 296 44 [bb 44] NOTE_INSN_BASIC_BLOCK)

(insn 296 295 297 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 135 [ save_eptr.165 ])) -1 (nil)
    (nil))

(insn 297 296 694 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 134 [ save_filt.166 ]))) -1 (nil)
    (nil))

(jump_insn 694 297 299 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 44, registers live:
 (nil)

(barrier 299 694 745)

(note/s 745 299 300 113 "" NOTE_INSN_DELETED_LABEL 216)

;; Start of basic block 45, registers live: (nil)
(code_label/s 300 745 301 45 170 "" [2 uses])

(note 301 300 302 45 [bb 45] NOTE_INSN_BASIC_BLOCK)

(insn 302 301 303 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 100 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 303 302 304 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 99 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 304 303 305 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 99 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 305 304 676 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 100 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 676 305 307 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 45, registers live:
 (nil)

(barrier 307 676 308)

;; Start of basic block 46, registers live: (nil)
(code_label 308 307 309 46 156 "" [2 uses])

(note 309 308 310 46 [bb 46] NOTE_INSN_BASIC_BLOCK)

(note 310 309 311 46 ("./CppUtilities.cc") 167)

(insn 311 310 312 46 ./CppUtilities.cc:167 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 140 [ D.36162 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 312 311 642 46 ./CppUtilities.cc:167 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 429)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 46, registers live:
 (nil)

;; Start of basic block 47, registers live: (nil)
(note 642 312 313 47 [bb 47] NOTE_INSN_BASIC_BLOCK)

(insn 313 642 314 47 ./CppUtilities.cc:167 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 139 [ ind1 ])
            (const_int -1 [0xffffffffffffffff]))) -1 (nil)
    (nil))

(jump_insn 314 313 316 47 ./CppUtilities.cc:167 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 429)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 47, registers live:
 (nil)

;; Start of basic block 48, registers live: (nil)
(note 316 314 317 48 [bb 48] NOTE_INSN_BASIC_BLOCK)

(note 317 316 318 48 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1908)

(insn 318 317 319 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 167)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 319 318 320 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 2 cx)
        (reg/v:DI 138 [ ind2 ])) -1 (nil)
    (nil))

(insn 320 319 321 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 1 dx)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 321 320 322 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 4 si)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(insn 322 321 323 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 5 di)
        (reg:DI 167)) -1 (nil)
    (nil))

(call_insn 323 322 325 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSsmm") [flags 0x41] <function_decl 0x2b4c42b2a700 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 48, registers live:
 (nil)

;; Start of basic block 49, registers live: (nil)
(note 325 323 326 49 [bb 49] NOTE_INSN_BASIC_BLOCK)

(note 326 325 327 49 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 327 326 328 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 168)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 328 327 329 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 168)) -1 (nil)
    (nil))

(insn 329 328 330 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 142 [ base.102 ])) -1 (nil)
    (nil))

(call_insn 330 329 643 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 49, registers live:
 (nil)

;; Start of basic block 50, registers live: (nil)
(note 643 330 331 50 [bb 50] NOTE_INSN_BASIC_BLOCK)

(insn 331 643 333 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 125 [ D.38154 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 333 331 334 50 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 334 333 335 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 87 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -48 [0xffffffffffffffd0])) [11 D.36150._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 335 334 336 50 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 336 335 337 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 88 [ D.40005 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 337 336 338 50 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 338 337 339 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 87 [ this ])
            (reg:DI 88 [ D.40005 ]))) -1 (nil)
    (nil))

(jump_insn 339 338 341 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 50, registers live:
 (nil)

;; Start of basic block 51, registers live: (nil)
(note 341 339 342 51 [bb 51] NOTE_INSN_BASIC_BLOCK)

(note 342 341 343 51 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 343 342 344 51 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 169)
                (plus:DI (reg:DI 87 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 344 343 345 51 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 345 344 346 51 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 169)) -1 (nil)
    (nil))

(call_insn 346 345 644 51 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 28 [0x1c])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 51, registers live:
 (nil)

;; Start of basic block 52, registers live: (nil)
(note 644 346 347 52 [bb 52] NOTE_INSN_BASIC_BLOCK)

(insn 347 644 349 52 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 89 [ D.40013 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 349 347 350 52 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 89 [ D.40013 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 350 349 352 52 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 8122 [0x1fba])
        (nil)))
;; End of basic block 52, registers live:
 (nil)

;; Start of basic block 53, registers live: (nil)
(note 352 350 353 53 [bb 53] NOTE_INSN_BASIC_BLOCK)

(note 353 352 354 53 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 354 353 355 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 170)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -5 [0xfffffffffffffffb])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 355 354 356 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 170)) -1 (nil)
    (nil))

(insn 356 355 357 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 87 [ this ])) -1 (nil)
    (nil))

(call_insn 357 356 358 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 358 357 359 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref:DI 623)) 548 {jump} (nil)
    (nil))
;; End of basic block 53, registers live:
 (nil)

(barrier 359 358 733)

;; Start of basic block 54, registers live: (nil)
(code_label/s 733 359 736 54 213 "" [1 uses])

(note 736 733 734 54 [bb 54] NOTE_INSN_BASIC_BLOCK)

(insn 734 736 735 54 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 735 734 360 54 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 360 735 362 54 "" NOTE_INSN_DELETED_LABEL 176)

(insn 362 360 363 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 91 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 363 362 364 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 90 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 364 363 365 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 90 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 365 364 737 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 91 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 737 365 368 54 "" NOTE_INSN_DELETED_LABEL 214)

(note/s 368 737 370 54 "" NOTE_INSN_DELETED_LABEL 177)

(insn 370 368 371 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 93 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 371 370 372 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 92 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 372 371 373 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 92 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 373 372 672 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 93 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 672 373 375 54 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 54, registers live:
 (nil)

(barrier 375 672 777)

;; Start of basic block 57, registers live: (nil)
(code_label/s 777 375 780 57 224 "" [1 uses])

(note 780 777 778 57 [bb 57] NOTE_INSN_BASIC_BLOCK)

(insn 778 780 779 57 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 779 778 376 57 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 376 779 378 57 "" NOTE_INSN_DELETED_LABEL 178)

(insn 378 376 379 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 132 [ save_filt.168 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 379 378 380 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 133 [ save_eptr.167 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(note 380 379 381 57 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 381 380 382 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 80 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -48 [0xffffffffffffffd0])) [11 D.36150._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 382 381 383 57 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 383 382 384 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 81 [ D.40061 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 384 383 385 57 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 385 384 386 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 80 [ this ])
            (reg:DI 81 [ D.40061 ]))) -1 (nil)
    (nil))

(jump_insn 386 385 388 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 415)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 57, registers live:
 (nil)

;; Start of basic block 59, registers live: (nil)
(note 388 386 389 59 [bb 59] NOTE_INSN_BASIC_BLOCK)

(note 389 388 390 59 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 390 389 391 59 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 171)
                (plus:DI (reg:DI 80 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 391 390 392 59 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 392 391 393 59 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 171)) -1 (nil)
    (nil))

(call_insn 393 392 645 59 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 32 [0x20])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 59, registers live:
 (nil)

;; Start of basic block 60, registers live: (nil)
(note 645 393 394 60 [bb 60] NOTE_INSN_BASIC_BLOCK)

(insn 394 645 396 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 82 [ D.40069 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 396 394 397 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 82 [ D.40069 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 397 396 399 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 415)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 60, registers live:
 (nil)

;; Start of basic block 61, registers live: (nil)
(note 399 397 400 61 [bb 61] NOTE_INSN_BASIC_BLOCK)

(note 400 399 401 61 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 401 400 402 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 172)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -4 [0xfffffffffffffffc])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 402 401 403 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 172)) -1 (nil)
    (nil))

(insn 403 402 404 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 80 [ this ])) -1 (nil)
    (nil))

(call_insn 404 403 405 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 405 404 406 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 415)) -1 (nil)
    (nil))
;; End of basic block 61, registers live:
 (nil)

(barrier 406 405 725)

;; Start of basic block 62, registers live: (nil)
(code_label/s 725 406 728 62 211 "" [1 uses])

(note 728 725 726 62 [bb 62] NOTE_INSN_BASIC_BLOCK)

(insn 726 728 727 62 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 727 726 407 62 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 407 727 409 62 "" NOTE_INSN_DELETED_LABEL 182)

(insn 409 407 410 62 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 84 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 410 409 411 62 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 83 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 411 410 412 62 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 83 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 412 411 666 62 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 84 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 666 412 414 62 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 421)) -1 (nil)
    (nil))
;; End of basic block 62, registers live:
 (nil)

(barrier 414 666 415)

;; Start of basic block 64, registers live: (nil)
(code_label 415 414 416 64 179 "" [3 uses])

(note 416 415 417 64 [bb 64] NOTE_INSN_BASIC_BLOCK)

(insn 417 416 418 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 133 [ save_eptr.167 ])) -1 (nil)
    (nil))

(insn 418 417 692 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 132 [ save_filt.168 ]))) -1 (nil)
    (nil))

(jump_insn 692 418 420 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 64, registers live:
 (nil)

(barrier 420 692 729)

(note/s 729 420 421 109 "" NOTE_INSN_DELETED_LABEL 212)

;; Start of basic block 65, registers live: (nil)
(code_label/s 421 729 422 65 183 "" [2 uses])

(note 422 421 423 65 [bb 65] NOTE_INSN_BASIC_BLOCK)

(insn 423 422 424 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 86 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 424 423 425 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 85 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 425 424 426 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 85 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 426 425 668 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 86 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 668 426 428 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 65, registers live:
 (nil)

(barrier 428 668 429)

;; Start of basic block 66, registers live: (nil)
(code_label 429 428 430 66 171 "" [2 uses])

(note 430 429 431 66 [bb 66] NOTE_INSN_BASIC_BLOCK)

(note 431 430 432 66 ("./CppUtilities.cc") 171)

(insn 432 431 433 66 ./CppUtilities.cc:171 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 141 [ D.36161 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 433 432 646 66 ./CppUtilities.cc:171 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 66, registers live:
 (nil)

;; Start of basic block 67, registers live: (nil)
(note 646 433 434 67 [bb 67] NOTE_INSN_BASIC_BLOCK)

(insn 434 646 435 67 ./CppUtilities.cc:171 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:DI 138 [ ind2 ])
            (const_int -1 [0xffffffffffffffff]))) -1 (nil)
    (nil))

(jump_insn 435 434 437 67 ./CppUtilities.cc:171 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 67, registers live:
 (nil)

;; Start of basic block 68, registers live: (nil)
(note 437 435 438 68 [bb 68] NOTE_INSN_BASIC_BLOCK)

(note 438 437 439 68 ("./CppUtilities.cc") 173)

(insn 439 438 440 68 ./CppUtilities.cc:173 (parallel [
            (set (reg/v:DI 79 [ __pos ])
                (plus:DI (reg/v:DI 139 [ ind1 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 440 439 441 68 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 305)

(insn 441 440 442 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (reg/f:DI 173)
        (mem/s/f:DI (reg/v/f:DI 144 [ file ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])) -1 (nil)
    (nil))

(insn 442 441 443 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (reg:CC 17 flags)
        (compare:CC (reg/v:DI 79 [ __pos ])
            (mem/s:DI (plus:DI (reg/f:DI 173)
                    (const_int -24 [0xffffffffffffffe8])) [19 <variable>.D.17834._M_length+0 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 443 442 445 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:305 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 450)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 68, registers live:
 (nil)

;; Start of basic block 69, registers live: (nil)
(note 445 443 446 69 [bb 69] NOTE_INSN_BASIC_BLOCK)

(note 446 445 447 69 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 306)

(insn 447 446 448 69 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:306 (set (reg:DI 5 di)
        (symbol_ref/f:DI ("*.LC7") [flags 0x2] <string_cst 0x2b4c43cf0380>)) -1 (nil)
    (nil))

(call_insn 448 447 449 69 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:306 (call (mem:QI (symbol_ref:DI ("_ZSt20__throw_out_of_rangePKc") [flags 0x41] <function_decl 0x2b4c423ce000 __throw_out_of_range>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 1 [0x1])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 69, registers live:
 (nil)

(barrier 449 448 450)

;; Start of basic block 70, registers live: (nil)
(code_label 450 449 451 70 185 "" [1 uses])

(note 451 450 452 70 [bb 70] NOTE_INSN_BASIC_BLOCK)

(note 452 451 453 70 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 1908)

(insn 453 452 454 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (parallel [
            (set (reg:DI 174)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -32 [0xffffffffffffffe0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 454 453 455 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 2 cx)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 455 454 456 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 1 dx)
        (reg/v:DI 79 [ __pos ])) -1 (nil)
    (nil))

(insn 456 455 457 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 4 si)
        (reg/v/f:DI 144 [ file ])) -1 (nil)
    (nil))

(insn 457 456 458 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (set (reg:DI 5 di)
        (reg:DI 174)) -1 (nil)
    (nil))

(call_insn 458 457 460 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:1908 (call (mem:QI (symbol_ref:DI ("_ZNSsC1ERKSsmm") [flags 0x41] <function_decl 0x2b4c42b2a700 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 70, registers live:
 (nil)

;; Start of basic block 71, registers live: (nil)
(note 460 458 461 71 [bb 71] NOTE_INSN_BASIC_BLOCK)

(note 461 460 462 71 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 462 461 463 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 175)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -32 [0xffffffffffffffe0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 463 462 464 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 175)) -1 (nil)
    (nil))

(insn 464 463 465 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 142 [ base.102 ])) -1 (nil)
    (nil))

(call_insn 465 464 647 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 8 [0x8])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 71, registers live:
 (nil)

;; Start of basic block 72, registers live: (nil)
(note 647 465 466 72 [bb 72] NOTE_INSN_BASIC_BLOCK)

(insn 466 647 468 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 124 [ D.38159 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 468 466 469 72 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 469 468 470 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 72 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -32 [0xffffffffffffffe0])) [11 D.36151._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 470 469 471 72 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 471 470 472 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 73 [ D.40139 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 472 471 473 72 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 473 472 474 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 72 [ this ])
            (reg:DI 73 [ D.40139 ]))) -1 (nil)
    (nil))

(jump_insn 474 473 476 72 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 72, registers live:
 (nil)

;; Start of basic block 73, registers live: (nil)
(note 476 474 477 73 [bb 73] NOTE_INSN_BASIC_BLOCK)

(note 477 476 478 73 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 478 477 479 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 176)
                (plus:DI (reg:DI 72 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 479 478 480 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 480 479 481 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 176)) -1 (nil)
    (nil))

(call_insn 481 480 648 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 36 [0x24])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 73, registers live:
 (nil)

;; Start of basic block 74, registers live: (nil)
(note 648 481 482 74 [bb 74] NOTE_INSN_BASIC_BLOCK)

(insn 482 648 484 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 74 [ D.40147 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 484 482 485 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 74 [ D.40147 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 485 484 487 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 623)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 8122 [0x1fba])
        (nil)))
;; End of basic block 74, registers live:
 (nil)

;; Start of basic block 75, registers live: (nil)
(note 487 485 488 75 [bb 75] NOTE_INSN_BASIC_BLOCK)

(note 488 487 489 75 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 489 488 490 75 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 177)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -3 [0xfffffffffffffffd])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 490 489 491 75 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 177)) -1 (nil)
    (nil))

(insn 491 490 492 75 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 72 [ this ])) -1 (nil)
    (nil))

(call_insn 492 491 493 75 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 493 492 494 75 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref:DI 623)) 548 {jump} (nil)
    (nil))
;; End of basic block 75, registers live:
 (nil)

(barrier 494 493 717)

;; Start of basic block 76, registers live: (nil)
(code_label/s 717 494 720 76 209 "" [1 uses])

(note 720 717 718 76 [bb 76] NOTE_INSN_BASIC_BLOCK)

(insn 718 720 719 76 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 719 718 495 76 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 495 719 497 76 "" NOTE_INSN_DELETED_LABEL 190)

(insn 497 495 498 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 76 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 498 497 499 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 75 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 499 498 500 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 75 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 500 499 721 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 76 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 721 500 503 76 "" NOTE_INSN_DELETED_LABEL 210)

(note/s 503 721 505 76 "" NOTE_INSN_DELETED_LABEL 191)

(insn 505 503 506 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 78 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 506 505 507 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 77 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 507 506 508 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 77 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 508 507 664 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 78 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 664 508 510 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 76, registers live:
 (nil)

(barrier 510 664 773)

;; Start of basic block 79, registers live: (nil)
(code_label/s 773 510 776 79 223 "" [1 uses])

(note 776 773 774 79 [bb 79] NOTE_INSN_BASIC_BLOCK)

(insn 774 776 775 79 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 775 774 511 79 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 511 775 513 79 "" NOTE_INSN_DELETED_LABEL 192)

(insn 513 511 514 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 130 [ save_filt.170 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 514 513 515 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 131 [ save_eptr.169 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(note 515 514 516 79 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 516 515 517 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 65 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -32 [0xffffffffffffffe0])) [11 D.36151._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 517 516 518 79 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 518 517 519 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 66 [ D.40195 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 519 518 520 79 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 520 519 521 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 65 [ this ])
            (reg:DI 66 [ D.40195 ]))) -1 (nil)
    (nil))

(jump_insn 521 520 523 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 550)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 79, registers live:
 (nil)

;; Start of basic block 81, registers live: (nil)
(note 523 521 524 81 [bb 81] NOTE_INSN_BASIC_BLOCK)

(note 524 523 525 81 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 525 524 526 81 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 178)
                (plus:DI (reg:DI 65 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 526 525 527 81 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 527 526 528 81 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 178)) -1 (nil)
    (nil))

(call_insn 528 527 649 81 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 40 [0x28])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 81, registers live:
 (nil)

;; Start of basic block 82, registers live: (nil)
(note 649 528 529 82 [bb 82] NOTE_INSN_BASIC_BLOCK)

(insn 529 649 531 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 67 [ D.40203 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 531 529 532 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 67 [ D.40203 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 532 531 534 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 550)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 82, registers live:
 (nil)

;; Start of basic block 83, registers live: (nil)
(note 534 532 535 83 [bb 83] NOTE_INSN_BASIC_BLOCK)

(note 535 534 536 83 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 536 535 537 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 179)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -2 [0xfffffffffffffffe])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 537 536 538 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 179)) -1 (nil)
    (nil))

(insn 538 537 539 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 65 [ this ])) -1 (nil)
    (nil))

(call_insn 539 538 540 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 540 539 541 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 550)) -1 (nil)
    (nil))
;; End of basic block 83, registers live:
 (nil)

(barrier 541 540 709)

;; Start of basic block 84, registers live: (nil)
(code_label/s 709 541 712 84 207 "" [1 uses])

(note 712 709 710 84 [bb 84] NOTE_INSN_BASIC_BLOCK)

(insn 710 712 711 84 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 711 710 542 84 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 542 711 544 84 "" NOTE_INSN_DELETED_LABEL 196)

(insn 544 542 545 84 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 69 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 545 544 546 84 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 68 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 546 545 547 84 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 68 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 547 546 658 84 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 69 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 658 547 549 84 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 556)) -1 (nil)
    (nil))
;; End of basic block 84, registers live:
 (nil)

(barrier 549 658 550)

;; Start of basic block 86, registers live: (nil)
(code_label 550 549 551 86 193 "" [3 uses])

(note 551 550 552 86 [bb 86] NOTE_INSN_BASIC_BLOCK)

(insn 552 551 553 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 131 [ save_eptr.169 ])) -1 (nil)
    (nil))

(insn 553 552 690 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 130 [ save_filt.170 ]))) -1 (nil)
    (nil))

(jump_insn 690 553 555 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 86, registers live:
 (nil)

(barrier 555 690 713)

(note/s 713 555 556 105 "" NOTE_INSN_DELETED_LABEL 208)

;; Start of basic block 87, registers live: (nil)
(code_label/s 556 713 557 87 197 "" [2 uses])

(note 557 556 558 87 [bb 87] NOTE_INSN_BASIC_BLOCK)

(insn 558 557 559 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 71 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 559 558 560 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 70 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 560 559 561 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 70 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 561 560 660 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 71 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 660 561 563 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 564)) -1 (nil)
    (nil))
;; End of basic block 87, registers live:
 (nil)

(barrier 563 660 789)

;; Start of basic block 88, registers live: (nil)
(code_label/s 789 563 792 88 227 "" [1 uses])

(note 792 789 790 88 [bb 88] NOTE_INSN_BASIC_BLOCK)

(insn 790 792 791 88 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 791 790 564 88 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 88, registers live:
 (nil)

;; Start of basic block 89, registers live: (nil)
(code_label/s 564 791 565 89 198 "" [13 uses])

(note 565 564 566 89 [bb 89] NOTE_INSN_BASIC_BLOCK)

(insn 566 565 567 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 128 [ save_filt.172 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 567 566 568 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 129 [ save_eptr.171 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(note 568 567 569 89 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 569 568 570 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 58 [ this ])
                (plus:DI (mem/s/f:DI (reg/f:DI 142 [ base.102 ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 570 569 571 89 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 571 570 572 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 59 [ D.40252 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 572 571 573 89 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 573 572 574 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 58 [ this ])
            (reg:DI 59 [ D.40252 ]))) -1 (nil)
    (nil))

(jump_insn 574 573 576 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 603)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 89, registers live:
 (nil)

;; Start of basic block 90, registers live: (nil)
(note 576 574 577 90 [bb 90] NOTE_INSN_BASIC_BLOCK)

(note 577 576 578 90 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 578 577 579 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 180)
                (plus:DI (reg:DI 58 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 579 578 580 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 580 579 581 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 180)) -1 (nil)
    (nil))

(call_insn 581 580 650 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 44 [0x2c])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 90, registers live:
 (nil)

;; Start of basic block 91, registers live: (nil)
(note 650 581 582 91 [bb 91] NOTE_INSN_BASIC_BLOCK)

(insn 582 650 584 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 60 [ D.40260 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 584 582 585 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 60 [ D.40260 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 585 584 587 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 603)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 91, registers live:
 (nil)

;; Start of basic block 92, registers live: (nil)
(note 587 585 588 92 [bb 92] NOTE_INSN_BASIC_BLOCK)

(note 588 587 589 92 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 589 588 590 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 181)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 590 589 591 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 181)) -1 (nil)
    (nil))

(insn 591 590 592 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 58 [ this ])) -1 (nil)
    (nil))

(call_insn 592 591 593 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 593 592 594 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 603)) -1 (nil)
    (nil))
;; End of basic block 92, registers live:
 (nil)

(barrier 594 593 701)

;; Start of basic block 93, registers live: (nil)
(code_label/s 701 594 704 93 205 "" [1 uses])

(note 704 701 702 93 [bb 93] NOTE_INSN_BASIC_BLOCK)

(insn 702 704 703 93 (set (reg:DI 153)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 703 702 595 93 (set (reg:DI 152)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 595 703 597 93 "" NOTE_INSN_DELETED_LABEL 202)

(insn 597 595 598 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 62 [ save_filt.160 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 598 597 599 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 61 [ save_eptr.159 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 599 598 600 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 61 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 600 599 653 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 62 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 653 600 602 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 609)) -1 (nil)
    (nil))
;; End of basic block 93, registers live:
 (nil)

(barrier 602 653 603)

;; Start of basic block 95, registers live: (nil)
(code_label 603 602 604 95 199 "" [3 uses])

(note 604 603 605 95 [bb 95] NOTE_INSN_BASIC_BLOCK)

(insn 605 604 606 95 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 129 [ save_eptr.171 ])) -1 (nil)
    (nil))

(insn 606 605 698 95 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 128 [ save_filt.172 ]))) -1 (nil)
    (nil))

(insn 698 606 699 95 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 153)) -1 (nil)
    (nil))

(call_insn 699 698 608 95 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 95, registers live:
 (nil)

(barrier 608 699 705)

(note/s 705 608 609 103 "" NOTE_INSN_DELETED_LABEL 206)

;; Start of basic block 96, registers live: (nil)
(code_label/s 609 705 610 96 203 "" [2 uses])

(note 610 609 611 96 [bb 96] NOTE_INSN_BASIC_BLOCK)

(insn 611 610 612 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 64 [ save_filt.162 ])
        (subreg:SI (reg:DI 152) 0)) -1 (nil)
    (nil))

(insn 612 611 613 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 63 [ save_eptr.161 ])
        (reg:DI 153)) -1 (nil)
    (nil))

(insn 613 612 614 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 153)
        (reg:DI 63 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 614 613 655 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 64 [ save_filt.162 ]))) -1 (nil)
    (nil))

(insn 655 614 656 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 153)) -1 (nil)
    (nil))

(call_insn 656 655 616 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 96, registers live:
 (nil)

(barrier 616 656 621)

(note 621 616 622 NOTE_INSN_FUNCTION_END)

(note 622 621 623 ("./CppUtilities.cc") 180)

;; Start of basic block 97, registers live: (nil)
(code_label 623 622 651 97 138 "" [14 uses])

(note 651 623 624 97 [bb 97] NOTE_INSN_BASIC_BLOCK)

(insn 624 651 625 97 ./CppUtilities.cc:180 (set (reg/i:DI 0 ax)
        (reg/f:DI 143 [ <result> ])) -1 (nil)
    (nil))

(insn 625 624 631 97 ./CppUtilities.cc:180 (set (reg/i:DI 0 ax)
        (reg/f:DI 143 [ <result> ])) -1 (nil)
    (nil))

(insn 631 625 0 97 ./CppUtilities.cc:180 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 97, registers live:
 (nil)


;; Function std::string baseFileName(const char*, bool) (_Z12baseFileNamePKcb)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 4 failed.
Forwarding edge 0->2 to 4 failed.
Forwarding edge 4->5 to 14 failed.
Deleted label in block 7.
Merged 8 and 9 without moving.
Deleted label in block 10.
Edge 14->20 redirected to 22
Deleted label in block 15.
Merged 16 and 17 without moving.
Edge 16->20 redirected to 22
Deleted label in block 18.
Redirecting jump 119 from 20 to 22.
Deleting block 20.
Merged 22 and 23 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 4 failed.
Forwarding edge 4->5 to 14 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 4 failed.
Merged 2 and 3 without moving.
Forwarding edge 4->5 to 15 failed.
Merged 6 and 7 without moving.
Merged 11 and 12 without moving.
Merged 19 and 20 without moving.
Merged 19 and 21 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->1 to 4 failed.
Forwarding edge 4->5 to 15 failed.
(note 1 0 9 ("./CppUtilities.cc") 139)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 ./CppUtilities.cc:139 (set (reg/f:DI 76 [ <result> ])
        (reg:DI 5 di [ D.40526 ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppUtilities.cc:139 (set (reg/v/f:DI 77 [ file ])
        (reg:DI 4 si [ file ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppUtilities.cc:139 (set (reg:SI 79)
        (reg:SI 1 dx [ withSuffix ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppUtilities.cc:139 (set (reg/v:QI 78 [ withSuffix ])
        (subreg:QI (reg:SI 79) 0)) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("./CppUtilities.cc") 141)

(insn 12 11 13 0 ./CppUtilities.cc:141 (parallel [
            (set (reg:DI 80)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -3 [0xfffffffffffffffd])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppUtilities.cc:141 (parallel [
            (set (reg:DI 81)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -16 [0xfffffffffffffff0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppUtilities.cc:141 (set (reg:DI 1 dx)
        (reg:DI 80)) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppUtilities.cc:141 (set (reg:DI 4 si)
        (reg/v/f:DI 77 [ file ])) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppUtilities.cc:141 (set (reg:DI 5 di)
        (reg:DI 81)) -1 (nil)
    (nil))

(call_insn 17 16 153 0 ./CppUtilities.cc:141 (call (mem:QI (symbol_ref:DI ("_ZNSsC1EPKcRKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2a100 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 153 17 18 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 18 153 19 1 ./CppUtilities.cc:141 (set (pc)
        (label_ref 28)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 19 18 195)

;; Start of basic block 2, registers live: (nil)
(code_label/s 195 19 198 2 256 "" [1 uses])

(note 198 195 196 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 196 198 197 2 (set (reg:DI 83)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 197 196 20 2 (set (reg:DI 82)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 20 197 22 2 "" NOTE_INSN_DELETED_LABEL 237)

(insn 22 20 23 2 ./CppUtilities.cc:141 (set (reg:SI 74 [ save_filt.174 ])
        (subreg:SI (reg:DI 82) 0)) -1 (nil)
    (nil))

(insn 23 22 24 2 ./CppUtilities.cc:141 (set (reg:DI 75 [ save_eptr.173 ])
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 24 23 25 2 ./CppUtilities.cc:141 (set (reg:DI 83)
        (reg:DI 75 [ save_eptr.173 ])) -1 (nil)
    (nil))

(insn 25 24 172 2 ./CppUtilities.cc:141 (set (reg:DI 82)
        (sign_extend:DI (reg:SI 74 [ save_filt.174 ]))) -1 (nil)
    (nil))

(insn 172 25 173 2 ./CppUtilities.cc:141 (set (reg:DI 5 di)
        (reg:DI 83)) -1 (nil)
    (nil))

(call_insn 173 172 27 2 ./CppUtilities.cc:141 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 2, registers live:
 (nil)

(barrier 27 173 28)

;; Start of basic block 4, registers live: (nil)
(code_label 28 27 29 4 236 "" [1 uses])

(note 29 28 30 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 30 29 31 4 ("./CppUtilities.cc") 143)

(insn 31 30 32 4 ./CppUtilities.cc:143 (parallel [
            (set (reg:SI 84)
                (zero_extend:SI (reg/v:QI 78 [ withSuffix ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 32 31 33 4 ./CppUtilities.cc:143 (parallel [
            (set (reg:DI 85)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -16 [0xfffffffffffffff0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 33 32 34 4 ./CppUtilities.cc:143 (set (reg:SI 1 dx)
        (reg:SI 84)) -1 (nil)
    (nil))

(insn 34 33 35 4 ./CppUtilities.cc:143 (set (reg:DI 4 si)
        (reg:DI 85)) -1 (nil)
    (nil))

(insn 35 34 36 4 ./CppUtilities.cc:143 (set (reg:DI 5 di)
        (reg/f:DI 76 [ <result> ])) -1 (nil)
    (nil))

(call_insn 36 35 154 4 ./CppUtilities.cc:143 (call (mem:QI (symbol_ref:DI ("_Z12baseFileNameRKSsb") [flags 0x3] <function_decl 0x2b4c4365c500 baseFileName>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 3 [0x3])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:SI 1 dx))
                (nil)))))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 154 36 37 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 37 154 38 5 ./CppUtilities.cc:143 (set (pc)
        (label_ref 92)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 38 37 191)

;; Start of basic block 6, registers live: (nil)
(code_label/s 191 38 194 6 255 "" [1 uses])

(note 194 191 192 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 192 194 193 6 (set (reg:DI 83)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 193 192 39 6 (set (reg:DI 82)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 39 193 41 6 "" NOTE_INSN_DELETED_LABEL 239)

(insn 41 39 42 6 ./CppUtilities.cc:143 (set (reg:SI 72 [ save_filt.176 ])
        (subreg:SI (reg:DI 82) 0)) -1 (nil)
    (nil))

(insn 42 41 43 6 ./CppUtilities.cc:143 (set (reg:DI 73 [ save_eptr.175 ])
        (reg:DI 83)) -1 (nil)
    (nil))

(note 43 42 44 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 44 43 45 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 65 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -16 [0xfffffffffffffff0])) [11 fileStr._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 45 44 46 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 46 45 47 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 66 [ D.40415 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 47 46 48 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 48 47 49 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 65 [ this ])
            (reg:DI 66 [ D.40415 ]))) -1 (nil)
    (nil))

(jump_insn 49 48 51 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 78)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 51 49 52 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 52 51 53 8 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 53 52 54 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 86)
                (plus:DI (reg:DI 65 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 54 53 55 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 55 54 56 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 86)) -1 (nil)
    (nil))

(call_insn 56 55 155 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 6 [0x6])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 155 56 57 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 57 155 59 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 67 [ D.40423 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 59 57 60 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 67 [ D.40423 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 60 59 62 9 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 78)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(note 62 60 63 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 63 62 64 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 64 63 65 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 87)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -2 [0xfffffffffffffffe])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 65 64 66 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 87)) -1 (nil)
    (nil))

(insn 66 65 67 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 65 [ this ])) -1 (nil)
    (nil))

(call_insn 67 66 68 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 68 67 69 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 78)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 69 68 183)

;; Start of basic block 11, registers live: (nil)
(code_label/s 183 69 186 11 253 "" [1 uses])

(note 186 183 184 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 184 186 185 11 (set (reg:DI 83)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 185 184 70 11 (set (reg:DI 82)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 70 185 72 11 "" NOTE_INSN_DELETED_LABEL 243)

(insn 72 70 73 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 69 [ save_filt.160 ])
        (subreg:SI (reg:DI 82) 0)) -1 (nil)
    (nil))

(insn 73 72 74 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 68 [ save_eptr.159 ])
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 74 73 75 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 83)
        (reg:DI 68 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 75 74 164 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 82)
        (sign_extend:DI (reg:SI 69 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 164 75 77 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 84)) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

(barrier 77 164 78)

;; Start of basic block 13, registers live: (nil)
(code_label 78 77 79 13 240 "" [3 uses])

(note 79 78 80 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 80 79 81 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 83)
        (reg:DI 73 [ save_eptr.175 ])) -1 (nil)
    (nil))

(insn 81 80 169 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 82)
        (sign_extend:DI (reg:SI 72 [ save_filt.176 ]))) -1 (nil)
    (nil))

(insn 169 81 170 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 83)) -1 (nil)
    (nil))

(call_insn 170 169 83 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 13, registers live:
 (nil)

(barrier 83 170 187)

(note/s 187 83 84 27 "" NOTE_INSN_DELETED_LABEL 254)

;; Start of basic block 14, registers live: (nil)
(code_label/s 84 187 85 14 244 "" [2 uses])

(note 85 84 86 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(insn 86 85 87 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 71 [ save_filt.162 ])
        (subreg:SI (reg:DI 82) 0)) -1 (nil)
    (nil))

(insn 87 86 88 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 70 [ save_eptr.161 ])
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 88 87 89 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 83)
        (reg:DI 70 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 89 88 166 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 82)
        (sign_extend:DI (reg:SI 71 [ save_filt.162 ]))) -1 (nil)
    (nil))

(insn 166 89 167 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 83)) -1 (nil)
    (nil))

(call_insn 167 166 91 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 14, registers live:
 (nil)

(barrier 91 167 92)

;; Start of basic block 15, registers live: (nil)
(code_label 92 91 93 15 238 "" [1 uses])

(note 93 92 94 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 94 93 95 15 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 95 94 96 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 58 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -16 [0xfffffffffffffff0])) [11 fileStr._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 96 95 97 15 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 97 96 98 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 59 [ D.40471 ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 98 97 99 15 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 99 98 100 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 58 [ this ])
            (reg:DI 59 [ D.40471 ]))) -1 (nil)
    (nil))

(jump_insn 100 99 102 15 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 143)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(note 102 100 103 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(note 103 102 104 16 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 104 103 105 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 88)
                (plus:DI (reg:DI 58 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 105 104 106 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 106 105 107 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 88)) -1 (nil)
    (nil))

(call_insn 107 106 156 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 10 [0xa])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(note 156 107 108 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(insn 108 156 110 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 60 [ D.40479 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 110 108 111 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 60 [ D.40479 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 111 110 113 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 143)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 8122 [0x1fba])
        (nil)))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(note 113 111 114 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(note 114 113 115 18 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 115 114 116 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 89)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 116 115 117 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 89)) -1 (nil)
    (nil))

(insn 117 116 118 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 58 [ this ])) -1 (nil)
    (nil))

(call_insn 118 117 119 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 119 118 120 18 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref:DI 143)) 548 {jump} (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

(barrier 120 119 175)

;; Start of basic block 19, registers live: (nil)
(code_label/s 175 120 178 19 251 "" [1 uses])

(note 178 175 176 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(insn 176 178 177 19 (set (reg:DI 83)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 177 176 121 19 (set (reg:DI 82)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 121 177 123 19 "" NOTE_INSN_DELETED_LABEL 248)

(insn 123 121 124 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 62 [ save_filt.160 ])
        (subreg:SI (reg:DI 82) 0)) -1 (nil)
    (nil))

(insn 124 123 125 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 61 [ save_eptr.159 ])
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 125 124 126 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 83)
        (reg:DI 61 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 126 125 179 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 82)
        (sign_extend:DI (reg:SI 62 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 179 126 133 19 "" NOTE_INSN_DELETED_LABEL 252)

(note/s 133 179 135 19 "" NOTE_INSN_DELETED_LABEL 249)

(insn 135 133 136 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 64 [ save_filt.162 ])
        (subreg:SI (reg:DI 82) 0)) -1 (nil)
    (nil))

(insn 136 135 137 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 63 [ save_eptr.161 ])
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 137 136 138 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 83)
        (reg:DI 63 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 138 137 161 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 82)
        (sign_extend:DI (reg:SI 64 [ save_filt.162 ]))) -1 (nil)
    (nil))

(insn 161 138 162 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 83)) -1 (nil)
    (nil))

(call_insn 162 161 140 19 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 19, registers live:
 (nil)

(barrier 140 162 141)

(note 141 140 142 NOTE_INSN_FUNCTION_END)

(note 142 141 143 ("./CppUtilities.cc") 144)

;; Start of basic block 22, registers live: (nil)
(code_label 143 142 157 22 235 "" [3 uses])

(note 157 143 144 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(insn 144 157 145 22 ./CppUtilities.cc:144 (set (reg/i:DI 0 ax)
        (reg/f:DI 76 [ <result> ])) -1 (nil)
    (nil))

(insn 145 144 151 22 ./CppUtilities.cc:144 (set (reg/i:DI 0 ax)
        (reg/f:DI 76 [ <result> ])) -1 (nil)
    (nil))

(insn 151 145 0 22 ./CppUtilities.cc:144 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 22, registers live:
 (nil)


;; Function std::string join(const char*, std::vector<int, std::allocator<int> >&) (_Z4joinPKcRSt6vectorIiSaIiEE)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Deleted label in block 5.
Merged 6 and 7 without moving.
Deleted label in block 8.
Merged 12 and 13 without moving.
Merged 15 and 16 without moving.
Merged 17 and 18 without moving.
Merged 21 and 22 without moving.
Deleted label in block 24.
Deleted label in block 27.
Deleted label in block 28.
Merged 30 and 31 without moving.
Deleted label in block 32.
Merged 33 and 34 without moving.
Deleted label in block 35.
Deleted label in block 39.
Merged 40 and 41 without moving.
Deleted label in block 42.
Merged 49 and 50 without moving.
Deleted label in block 51.
Merged 52 and 53 without moving.
Deleted label in block 54.
Deleted label in block 58.
Merged 59 and 60 without moving.
Deleted label in block 61.
Deleted label in block 69.
Merged 70 and 71 without moving.
Deleted label in block 72.
Deleted label in block 77.
Merged 78 and 79 without moving.
Deleted label in block 80.
Merged 84 and 85 without moving.
Merged 87 and 88 without moving.
Merged 90 and 91 without moving.
Deleted label in block 94.
Merged 95 and 96 without moving.
Deleted label in block 97.
Merged 101 and 102 without moving.
Merged 104 and 105 without moving.
Redirecting jump 806 from 109 to 110.
Merged 107 and 108 without moving.
Deleting block 109.
Merged 110 and 111 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

Merged 3 and 4 without moving.
Merged 8 and 9 without moving.
Merged 11 and 12 without moving.
Merged 14 and 15 without moving.
Merged 30 and 31 without moving.
Merged 30 and 32 without moving.
Merged 33 and 34 without moving.
Merged 38 and 39 without moving.
Merged 50 and 51 without moving.
Merged 50 and 52 without moving.
Merged 53 and 54 without moving.
Merged 58 and 59 without moving.
Merged 71 and 72 without moving.
Merged 80 and 81 without moving.
Merged 83 and 84 without moving.
Merged 94 and 95 without moving.
Merged 97 and 98 without moving.


try_optimize_cfg iteration 2

(note 1 0 12 ("./CppUtilities.cc") 40)

;; Start of basic block 0, registers live: (nil)
(note 12 1 7 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 7 12 8 0 ./CppUtilities.cc:40 (set (reg/f:DI 180 [ <result> ])
        (reg:DI 5 di [ D.41885 ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppUtilities.cc:40 (set (reg/v/f:DI 181 [ tok ])
        (reg:DI 4 si [ tok ])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppUtilities.cc:40 (set (reg/v/f:DI 182 [ v ])
        (reg:DI 1 dx [ v ])) -1 (nil)
    (nil))

(note 10 9 14 0 NOTE_INSN_FUNCTION_BEG)

(insn 14 10 15 0 ./CppUtilities.cc:40 (set (reg/f:DI 172 [ this ])
        (reg/v/f:DI 182 [ v ])) -1 (nil)
    (nil))

(note 15 14 16 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 809)

(insn 16 15 17 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 184)
        (mem/s/f:DI (reg/f:DI 172 [ this ]) [59 <variable>.D.35879._M_impl._M_start+0 S8 A64])) -1 (nil)
    (nil))

(insn 17 16 18 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 185)
        (mem/s/f:DI (plus:DI (reg/f:DI 172 [ this ])
                (const_int 8 [0x8])) [59 <variable>.D.35879._M_impl._M_finish+0 S8 A64])) -1 (nil)
    (nil))

(insn 18 17 19 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 183)
                (minus:DI (reg:DI 185)
                    (reg:DI 184)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (minus:DI (mem/s/f:DI (plus:DI (reg/f:DI 172 [ this ])
                    (const_int 8 [0x8])) [59 <variable>.D.35879._M_impl._M_finish+0 S8 A64])
            (mem/s/f:DI (reg/f:DI 172 [ this ]) [59 <variable>.D.35879._M_impl._M_start+0 S8 A64]))
        (nil)))

(insn 19 18 20 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (parallel [
            (set (reg:DI 186)
                (ashiftrt:DI (reg:DI 183)
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:DI (reg:DI 183)
            (const_int 4 [0x4]))
        (nil)))

(insn 20 19 21 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:809 (set (reg:DI 171 [ D.40586 ])
        (reg:DI 186)) -1 (nil)
    (nil))

(note 21 20 22 0 ("./CppUtilities.cc") 42)

(insn 22 21 23 0 ./CppUtilities.cc:42 (parallel [
            (set (reg/v:SI 176 [ n ])
                (plus:SI (subreg:SI (reg:DI 171 [ D.40586 ]) 0)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 23 22 24 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 413)

(insn 24 23 25 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (parallel [
            (set (reg:DI 187)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 25 24 26 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (parallel [
            (set (reg:DI 158 [ D.40600 ])
                (plus:DI (reg:DI 187)
                    (const_int 88 [0x58])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 26 25 27 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (set (reg:DI 157 [ this ])
        (reg:DI 158 [ D.40600 ])) -1 (nil)
    (nil))

(note 27 26 28 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h") 446)

(insn 28 27 29 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (reg:DI 5 di)
        (reg:DI 157 [ this ])) -1 (nil)
    (nil))

(call_insn 29 28 30 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_baseC2Ev") [flags 0x41] <function_decl 0x2b4c42c93600 __base_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 30 29 31 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:DI (reg:DI 157 [ this ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9basic_iosIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c430446e0 _ZTVSt9basic_iosIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 31 30 32 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:DI (plus:DI (reg:DI 157 [ this ])
                (const_int 216 [0xd8])) [46 <variable>._M_tie+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 32 31 33 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:QI (plus:DI (reg:DI 157 [ this ])
                (const_int 224 [0xe0])) [0 <variable>._M_fill+0 S1 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 33 32 34 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:QI (plus:DI (reg:DI 157 [ this ])
                (const_int 225 [0xe1])) [47 <variable>._M_fill_init+0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 34 33 35 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:DI (plus:DI (reg:DI 157 [ this ])
                (const_int 232 [0xe8])) [48 <variable>._M_streambuf+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 35 34 36 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:DI (plus:DI (reg:DI 157 [ this ])
                (const_int 240 [0xf0])) [49 <variable>._M_ctype+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 36 35 37 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:DI (plus:DI (reg:DI 157 [ this ])
                (const_int 248 [0xf8])) [50 <variable>._M_num_put+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 37 36 38 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:446 (set (mem/s:DI (plus:DI (reg:DI 157 [ this ])
                (const_int 256 [0x100])) [51 <variable>._M_num_get+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 38 37 39 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream") 341)

(insn 39 38 40 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (reg/f:DI 188)
        (const:DI (plus:DI (symbol_ref:DI ("_ZTTSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436d04d0 _ZTTSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 8 [0x8])))) -1 (nil)
    (nil))

(insn 40 39 41 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (reg:DI 169 [ D.40609 ])
        (mem/s/u/f:DI (reg/f:DI 188) [3 _ZTTSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE+8 S8 A64])) -1 (nil)
    (nil))

(insn 41 40 42 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -400 [0xfffffffffffffe70])) [6 s.D.34935._vptr.basic_ostream+0 S8 A64])
        (reg:DI 169 [ D.40609 ])) -1 (nil)
    (nil))

(insn 42 41 43 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (parallel [
            (set (reg:DI 170 [ this.129 ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 43 42 44 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (parallel [
            (set (reg:DI 168 [ D.40611 ])
                (plus:DI (reg:DI 169 [ D.40609 ])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 44 43 45 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (reg/f:DI 189)
        (const:DI (plus:DI (symbol_ref:DI ("_ZTTSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436d04d0 _ZTTSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 45 44 46 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (reg:DI 167 [ iftmp.131 ])
        (mem/s/u/f:DI (reg/f:DI 189) [3 _ZTTSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE+16 S8 A64])) -1 (nil)
    (nil))

(insn 46 45 47 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (reg:DI 190)
        (mem:DI (reg:DI 168 [ D.40611 ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 47 46 48 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:341 (set (mem/s:DI (plus:DI (reg:DI 170 [ this.129 ])
                (reg:DI 190)) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (reg:DI 167 [ iftmp.131 ])) -1 (nil)
    (nil))

(note 48 47 49 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 413)

(insn 49 48 50 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -400 [0xfffffffffffffe70])) [6 s.D.34935._vptr.basic_ostream+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436cad10 _ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 24 [0x18])))) -1 (nil)
    (nil))

(insn 50 49 51 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (parallel [
            (set (reg:DI 159 [ this.108 ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 51 50 52 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (parallel [
            (set (reg:DI 160 [ D.40620 ])
                (plus:DI (reg:DI 159 [ this.108 ])
                    (const_int 88 [0x58])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 52 51 53 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:413 (set (mem/s:DI (reg:DI 160 [ D.40620 ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436cad10 _ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 64 [0x40])))) -1 (nil)
    (nil))

(note 53 52 54 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 449)

(insn 54 53 55 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -392 [0xfffffffffffffe78])) [6 s._M_stringbuf.D.34785._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 55 54 56 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -384 [0xfffffffffffffe80])) [11 s._M_stringbuf.D.34785._M_in_beg+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 56 55 57 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -376 [0xfffffffffffffe88])) [11 s._M_stringbuf.D.34785._M_in_cur+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 57 56 58 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -368 [0xfffffffffffffe90])) [11 s._M_stringbuf.D.34785._M_in_end+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 58 57 59 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -360 [0xfffffffffffffe98])) [11 s._M_stringbuf.D.34785._M_out_beg+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 59 58 60 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -352 [0xfffffffffffffea0])) [11 s._M_stringbuf.D.34785._M_out_cur+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 60 59 61 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -344 [0xfffffffffffffea8])) [11 s._M_stringbuf.D.34785._M_out_end+0 S8 A64])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 61 60 62 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (parallel [
            (set (reg:DI 191)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 62 61 63 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (parallel [
            (set (reg:DI 192)
                (plus:DI (reg:DI 191)
                    (const_int 64 [0x40])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 63 62 64 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (set (reg:DI 5 di)
        (reg:DI 192)) -1 (nil)
    (nil))

(call_insn 64 63 65 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:449 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeC1Ev") [flags 0x41] <function_decl 0x2b4c42c1d100 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 65 64 66 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 101)

(insn 66 65 67 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:101 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -392 [0xfffffffffffffe78])) [6 s._M_stringbuf.D.34785._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436ba370 _ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 67 66 68 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:101 (set (mem/s/c:SI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -328 [0xfffffffffffffeb8])) [64 s._M_stringbuf._M_mode+0 S4 A64])
        (const_int 16 [0x10])) -1 (nil)
    (nil))

(note 68 67 69 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 182)

(insn 69 68 70 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:182 (set (reg:DI 156 [ this ])
        (symbol_ref:DI ("_ZNSs4_Rep20_S_empty_rep_storageE") [flags 0x40] <var_decl 0x2b4c42b2f370 _S_empty_rep_storage>)) -1 (nil)
    (nil))

(note 70 69 71 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 213)

(insn 71 70 72 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:213 (parallel [
            (set (reg/v/f:DI 155 [ __dat ])
                (plus:DI (reg:DI 156 [ this ])
                    (const_int 24 [0x18])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 72 71 73 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 257)

(insn 73 72 74 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:257 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -320 [0xfffffffffffffec0])) [11 s._M_stringbuf._M_string._M_dataplus._M_p+0 S8 A64])
        (reg/v/f:DI 155 [ __dat ])) -1 (nil)
    (nil))

(note 74 73 75 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 414)

(insn 75 74 76 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:414 (parallel [
            (set (reg:DI 193)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 76 75 77 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:414 (parallel [
            (set (reg:DI 194)
                (plus:DI (reg:DI 193)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 77 76 78 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:414 (set (reg:DI 4 si)
        (reg:DI 194)) -1 (nil)
    (nil))

(insn 78 77 79 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:414 (set (reg:DI 5 di)
        (reg:DI 160 [ D.40620 ])) -1 (nil)
    (nil))

(call_insn 79 78 81 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:414 (call (mem:QI (symbol_ref:DI ("_ZNSt9basic_iosIcSt11char_traitsIcEE4initEPSt15basic_streambufIcS1_E") [flags 0x41] <function_decl 0x2b4c43041e00 init>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 5 [0x5])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 81 79 82 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 82 81 83 1 ("./CppUtilities.cc") 46)

(insn 83 82 84 1 ./CppUtilities.cc:46 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 176 [ n ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 84 83 86 1 ./CppUtilities.cc:46 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 204)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 86 84 87 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 87 86 88 2 ./CppUtilities.cc:46 (set (reg:DI 63 [ prephitmp.870 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 88 87 89 2 ./CppUtilities.cc:46 (set (pc)
        (label_ref 234)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 89 88 1032)

;; Start of basic block 3, registers live: (nil)
(code_label/s 1032 89 1035 3 372 "" [1 uses])

(note 1035 1032 1033 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 1033 1035 1034 3 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1034 1033 90 3 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 90 1034 92 3 "" NOTE_INSN_DELETED_LABEL 266)

(insn 92 90 93 3 ./CppUtilities.cc:46 (set (reg:SI 162 [ save_filt.180 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 93 92 94 3 ./CppUtilities.cc:46 (set (reg:DI 161 [ save_eptr.179 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 94 93 95 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 95 94 96 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -392 [0xfffffffffffffe78])) [6 s._M_stringbuf.D.34785._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436ba370 _ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 96 95 97 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 97 96 98 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 147 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -320 [0xfffffffffffffec0])) [11 s._M_stringbuf._M_string._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 98 97 99 3 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 99 98 100 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 147 [ this ]))) -1 (nil)
    (nil))

(jump_insn 100 99 102 3 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 129)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 102 100 103 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 103 102 104 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 104 103 105 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 197)
                (plus:DI (reg:DI 147 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 105 104 106 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 106 105 107 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 197)) -1 (nil)
    (nil))

(call_insn 107 106 838 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 16 [0x10])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(note 838 107 108 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 108 838 110 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 148 [ D.40767 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 110 108 111 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 148 [ D.40767 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 111 110 113 6 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 129)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 113 111 114 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 114 113 115 7 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 115 114 116 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 198)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -10 [0xfffffffffffffff6])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 116 115 117 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 198)) -1 (nil)
    (nil))

(insn 117 116 118 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 147 [ this ])) -1 (nil)
    (nil))

(call_insn 118 117 119 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 119 118 120 7 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 129)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

(barrier 120 119 1020)

;; Start of basic block 8, registers live: (nil)
(code_label/s 1020 120 1023 8 369 "" [1 uses])

(note 1023 1020 1021 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(insn 1021 1023 1022 8 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1022 1021 121 8 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 121 1022 123 8 "" NOTE_INSN_DELETED_LABEL 270)

(insn 123 121 124 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 150 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 124 123 125 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 149 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 125 124 126 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 149 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 126 125 904 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 150 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 904 126 128 8 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 141)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

(barrier 128 904 129)

;; Start of basic block 10, registers live: (nil)
(code_label 129 128 130 10 267 "" [3 uses])

(note 130 129 131 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 131 130 132 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 132 131 133 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 199)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 133 132 134 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 145 [ this ])
                (plus:DI (reg:DI 199)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 134 133 135 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 200)

(insn 135 134 136 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (mem/s:DI (reg:DI 145 [ this ]) [6 <variable>._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 136 135 137 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (parallel [
            (set (reg:DI 200)
                (plus:DI (reg:DI 145 [ this ])
                    (const_int 56 [0x38])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 137 136 138 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 5 di)
        (reg:DI 200)) -1 (nil)
    (nil))

(call_insn 138 137 139 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeD1Ev") [flags 0x41] <function_decl 0x2b4c42c1d300 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 139 138 140 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (pc)
        (label_ref 168)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 140 139 1024)

(note/s 1024 140 141 138 "" NOTE_INSN_DELETED_LABEL 370)

;; Start of basic block 11, registers live: (nil)
(code_label/s 141 1024 142 11 272 "" [2 uses])

(note 142 141 143 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 143 142 144 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:SI 152 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 144 143 145 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 151 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 145 144 146 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 196)
        (reg:DI 151 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 146 145 1028 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 152 [ save_filt.162 ]))) -1 (nil)
    (nil))

(note/s 1028 146 149 11 "" NOTE_INSN_DELETED_LABEL 371)

(note/s 149 1028 151 11 "" NOTE_INSN_DELETED_LABEL 273)

(insn 151 149 152 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:SI 154 [ save_filt.190 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 152 151 153 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 153 [ save_eptr.189 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 153 152 154 11 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 154 153 155 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 201)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 155 154 156 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 146 [ this ])
                (plus:DI (reg:DI 201)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 156 155 157 11 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 200)

(insn 157 156 158 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (mem/s:DI (reg:DI 146 [ this ]) [6 <variable>._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 158 157 159 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (parallel [
            (set (reg:DI 202)
                (plus:DI (reg:DI 146 [ this ])
                    (const_int 56 [0x38])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 159 158 160 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 5 di)
        (reg:DI 202)) -1 (nil)
    (nil))

(call_insn 160 159 163 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeD1Ev") [flags 0x41] <function_decl 0x2b4c42c1d300 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 163 160 164 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 196)
        (reg:DI 153 [ save_eptr.189 ])) -1 (nil)
    (nil))

(insn 164 163 908 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 154 [ save_filt.190 ]))) -1 (nil)
    (nil))

(jump_insn 908 164 166 11 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (pc)
        (label_ref 174)) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

(barrier 166 908 167)

(note 167 166 168 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

;; Start of basic block 13, registers live: (nil)
(code_label 168 167 169 13 271 "" [1 uses])

(note 169 168 170 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 170 169 171 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (reg:DI 196)
        (reg:DI 161 [ save_eptr.179 ])) -1 (nil)
    (nil))

(insn 171 170 1036 13 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 162 [ save_filt.180 ]))) -1 (nil)
    (nil))
;; End of basic block 13, registers live:
 (nil)

(note/s 1036 171 174 141 "" NOTE_INSN_DELETED_LABEL 373)

;; Start of basic block 14, registers live: (nil)
(code_label/s 174 1036 175 14 275 "" [2 uses])

(note 175 174 176 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(insn 176 175 177 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (reg:SI 164 [ save_filt.182 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 177 176 178 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (reg:DI 163 [ save_eptr.181 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 178 177 179 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 414)

(insn 179 178 180 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:414 (parallel [
            (set (reg:DI 144 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 180 179 181 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream") 113)

(insn 181 180 182 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (reg:DI 144 [ this ]) [6 <variable>._vptr.basic_ostream+0 S8 A64])
        (reg:DI 169 [ D.40609 ])) -1 (nil)
    (nil))

(insn 182 181 183 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 203)
        (mem:DI (reg:DI 168 [ D.40611 ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 183 182 186 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (plus:DI (reg:DI 144 [ this ])
                (reg:DI 203)) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (reg:DI 167 [ iftmp.131 ])) -1 (nil)
    (nil))

(insn 186 183 187 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 196)
        (reg:DI 163 [ save_eptr.181 ])) -1 (nil)
    (nil))

(insn 187 186 1040 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 164 [ save_filt.182 ]))) -1 (nil)
    (nil))

(note/s 1040 187 190 14 "" NOTE_INSN_DELETED_LABEL 374)

(note/s 190 1040 192 14 "" NOTE_INSN_DELETED_LABEL 277)

(insn 192 190 193 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:SI 166 [ save_filt.184 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 193 192 194 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 165 [ save_eptr.183 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 194 193 195 14 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h") 268)

(insn 195 194 196 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (mem/s:DI (reg:DI 157 [ this ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9basic_iosIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c430446e0 _ZTVSt9basic_iosIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 196 195 197 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 157 [ this ])) -1 (nil)
    (nil))

(call_insn 197 196 200 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_baseD2Ev") [flags 0x41] <function_decl 0x2b4c42c93900 __base_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 8 [0x8])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 200 197 201 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 196)
        (reg:DI 165 [ save_eptr.183 ])) -1 (nil)
    (nil))

(insn 201 200 914 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 166 [ save_filt.184 ]))) -1 (nil)
    (nil))

(insn 914 201 915 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 196)) -1 (nil)
    (nil))

(call_insn 915 914 203 14 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 14, registers live:
 (nil)

(barrier 203 915 204)

;; Start of basic block 16, registers live: (nil)
(code_label 204 203 205 16 263 "" [1 uses])

(note 205 204 206 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(insn 206 205 207 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg/v:SI 175 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 207 206 208 16 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 62 [ ivtmp.874 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(code_label 208 207 209 17 279 "" [1 uses])

(note 209 208 210 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(note 210 209 211 17 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 704)

(insn 211 210 212 17 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:704 (parallel [
            (set (reg:DI 143 [ D.40848 ])
                (plus:DI (reg:DI 62 [ ivtmp.874 ])
                    (mem/s/f:DI (reg/v/f:DI 182 [ v ]) [59 <variable>.D.35879._M_impl._M_start+0 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 212 211 213 17 ("./CppUtilities.cc") 47)

(insn 213 212 214 17 ./CppUtilities.cc:47 (set (reg:SI 179 [ D.35951 ])
        (mem:SI (reg:DI 143 [ D.40848 ]) [16 S4 A32])) -1 (nil)
    (nil))

(insn 214 213 215 17 ./CppUtilities.cc:47 (parallel [
            (set (reg:DI 204)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 215 214 216 17 ./CppUtilities.cc:47 (set (reg:SI 4 si)
        (reg:SI 179 [ D.35951 ])) -1 (nil)
    (nil))

(insn 216 215 217 17 ./CppUtilities.cc:47 (set (reg:DI 5 di)
        (reg:DI 204)) -1 (nil)
    (nil))

(call_insn 217 216 839 17 ./CppUtilities.cc:47 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEi") [flags 0x41] <function_decl 0x2b4c43357600 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(note 839 217 218 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(insn 218 839 220 18 ./CppUtilities.cc:47 (set (reg:DI 178 [ D.35952 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 220 218 221 18 ./CppUtilities.cc:47 (set (reg:DI 4 si)
        (reg/v/f:DI 181 [ tok ])) -1 (nil)
    (nil))

(insn 221 220 222 18 ./CppUtilities.cc:47 (set (reg:DI 5 di)
        (reg:DI 178 [ D.35952 ])) -1 (nil)
    (nil))

(call_insn 222 221 223 18 ./CppUtilities.cc:47 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b4c43370400 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 18, registers live:
 (nil)

;; Start of basic block 19, registers live: (nil)
(note 223 222 224 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(note 224 223 225 19 ("./CppUtilities.cc") 46)

(insn 225 224 226 19 ./CppUtilities.cc:46 (parallel [
            (set (reg/v:SI 175 [ i ])
                (plus:SI (reg/v:SI 175 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 226 225 227 19 ./CppUtilities.cc:46 (parallel [
            (set (reg:DI 62 [ ivtmp.874 ])
                (plus:DI (reg:DI 62 [ ivtmp.874 ])
                    (const_int 4 [0x4])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 227 226 228 19 ./CppUtilities.cc:46 (parallel [
            (set (reg:SI 205)
                (plus:SI (subreg:SI (reg:DI 171 [ D.40586 ]) 0)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 228 227 229 19 ./CppUtilities.cc:46 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 175 [ i ])
            (reg:SI 205))) -1 (nil)
    (nil))

(jump_insn 229 228 231 19 ./CppUtilities.cc:46 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 208)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 19, registers live:
 (nil)

;; Start of basic block 20, registers live: (nil)
(note 231 229 232 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(insn 232 231 233 20 ./CppUtilities.cc:46 (set (reg:DI 206)
        (sign_extend:DI (reg/v:SI 175 [ i ]))) -1 (nil)
    (nil))

(insn 233 232 234 20 ./CppUtilities.cc:46 (parallel [
            (set (reg:DI 63 [ prephitmp.870 ])
                (ashift:DI (reg:DI 206)
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 20, registers live:
 (nil)

;; Start of basic block 21, registers live: (nil)
(code_label 234 233 235 21 265 "" [1 uses])

(note 235 234 236 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(note 236 235 237 21 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h") 704)

(insn 237 236 238 21 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/stl_iterator.h:704 (parallel [
            (set (reg:DI 142 [ D.40879 ])
                (plus:DI (reg:DI 63 [ prephitmp.870 ])
                    (mem/s/f:DI (reg/v/f:DI 182 [ v ]) [59 <variable>.D.35879._M_impl._M_start+0 S8 A64])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 238 237 239 21 ("./CppUtilities.cc") 48)

(insn 239 238 240 21 ./CppUtilities.cc:48 (set (reg:SI 177 [ D.35954 ])
        (mem:SI (reg:DI 142 [ D.40879 ]) [16 S4 A32])) -1 (nil)
    (nil))

(insn 240 239 241 21 ./CppUtilities.cc:48 (parallel [
            (set (reg:DI 207)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 241 240 242 21 ./CppUtilities.cc:48 (set (reg:SI 4 si)
        (reg:SI 177 [ D.35954 ])) -1 (nil)
    (nil))

(insn 242 241 243 21 ./CppUtilities.cc:48 (set (reg:DI 5 di)
        (reg:DI 207)) -1 (nil)
    (nil))

(call_insn 243 242 244 21 ./CppUtilities.cc:48 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEi") [flags 0x41] <function_decl 0x2b4c43357600 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 1 [0x1])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 21, registers live:
 (nil)

;; Start of basic block 22, registers live: (nil)
(note 244 243 245 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(note 245 244 246 22 ("./CppUtilities.cc") 50)

(insn 246 245 247 22 ./CppUtilities.cc:50 (set (reg/f:DI 98 [ this ])
        (reg/f:DI 180 [ <result> ])) -1 (nil)
    (nil))

(note 247 246 248 22 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 257)

(insn 248 247 249 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:257 (set (mem/s:DI (reg/f:DI 98 [ this ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])
        (reg/v/f:DI 155 [ __dat ])) -1 (nil)
    (nil))

(note 249 248 250 22 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 515)

(insn 250 249 251 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:515 (set (reg:DI 127 [ D.40954 ])
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -352 [0xfffffffffffffea0])) [11 s._M_stringbuf.D.34785._M_out_cur+0 S8 A64])) -1 (nil)
    (nil))

(note 251 250 252 22 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 131)

(insn 252 251 253 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:131 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 127 [ D.40954 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 253 252 255 22 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:131 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 499)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 913 [0x391])
        (nil)))
;; End of basic block 22, registers live:
 (nil)

;; Start of basic block 23, registers live: (nil)
(note 255 253 256 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(note 256 255 257 23 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 471)

(insn 257 256 258 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:471 (set (reg:DI 128 [ D.40961 ])
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -368 [0xfffffffffffffe90])) [11 s._M_stringbuf.D.34785._M_in_end+0 S8 A64])) -1 (nil)
    (nil))

(note 258 257 259 23 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 134)

(insn 259 258 260 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:134 (set (reg:CC 17 flags)
        (compare:CC (reg:DI 127 [ D.40954 ])
            (reg:DI 128 [ D.40961 ]))) -1 (nil)
    (nil))

(jump_insn 260 259 262 23 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:134 (set (pc)
        (if_then_else (leu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 380)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 23, registers live:
 (nil)

;; Start of basic block 24, registers live: (nil)
(note 262 260 263 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(note 263 262 264 24 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 135)

(insn 264 263 265 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (parallel [
            (set (reg:DI 208)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -8 [0xfffffffffffffff8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 265 264 266 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (set (reg:DI 209)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -360 [0xfffffffffffffe98])) [11 s._M_stringbuf.D.34785._M_out_beg+0 S8 A64])) -1 (nil)
    (nil))

(insn 266 265 267 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (parallel [
            (set (reg:DI 210)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -32 [0xffffffffffffffe0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 267 266 268 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (set (reg:DI 2 cx)
        (reg:DI 208)) -1 (nil)
    (nil))

(insn 268 267 269 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (set (reg:DI 1 dx)
        (reg:DI 127 [ D.40954 ])) -1 (nil)
    (nil))

(insn 269 268 270 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (set (reg:DI 4 si)
        (reg:DI 209)) -1 (nil)
    (nil))

(insn 270 269 271 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (set (reg:DI 5 di)
        (reg:DI 210)) -1 (nil)
    (nil))

(call_insn 271 270 272 24 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:135 (call (mem:QI (symbol_ref/i:DI ("_ZNSsC1IPcEET_S1_RKSaIcE") [flags 0x1] <function_decl 0x2b4c43819d00 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 20 [0x14])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(note 272 271 273 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(note 273 272 274 25 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 274 273 275 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 211)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -32 [0xffffffffffffffe0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 275 274 276 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 211)) -1 (nil)
    (nil))

(insn 276 275 277 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 98 [ this ])) -1 (nil)
    (nil))

(call_insn 277 276 840 25 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 21 [0x15])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 25, registers live:
 (nil)

;; Start of basic block 26, registers live: (nil)
(note 840 277 278 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(insn 278 840 280 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 139 [ D.40973 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 280 278 281 26 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 281 280 282 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 121 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -32 [0xffffffffffffffe0])) [11 D.40971._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 282 281 283 26 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 283 282 284 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 121 [ this ]))) -1 (nil)
    (nil))

(jump_insn 284 283 286 26 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 696)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 26, registers live:
 (nil)

;; Start of basic block 27, registers live: (nil)
(note 286 284 287 27 [bb 27] NOTE_INSN_BASIC_BLOCK)

(note 287 286 288 27 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 288 287 289 27 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 212)
                (plus:DI (reg:DI 121 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 289 288 290 27 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 290 289 291 27 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 212)) -1 (nil)
    (nil))

(call_insn 291 290 841 27 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 32 [0x20])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 27, registers live:
 (nil)

;; Start of basic block 28, registers live: (nil)
(note 841 291 292 28 [bb 28] NOTE_INSN_BASIC_BLOCK)

(insn 292 841 294 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 122 [ D.41059 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 294 292 295 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 122 [ D.41059 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 295 294 297 28 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 696)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 28, registers live:
 (nil)

;; Start of basic block 29, registers live: (nil)
(note 297 295 298 29 [bb 29] NOTE_INSN_BASIC_BLOCK)

(note 298 297 299 29 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 299 298 300 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 213)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -7 [0xfffffffffffffff9])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 300 299 301 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 213)) -1 (nil)
    (nil))

(insn 301 300 302 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 121 [ this ])) -1 (nil)
    (nil))

(call_insn 302 301 303 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 303 302 304 29 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 696)) -1 (nil)
    (nil))
;; End of basic block 29, registers live:
 (nil)

(barrier 304 303 992)

;; Start of basic block 30, registers live: (nil)
(code_label/s 992 304 995 30 362 "" [1 uses])

(note 995 992 993 30 [bb 30] NOTE_INSN_BASIC_BLOCK)

(insn 993 995 994 30 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 994 993 305 30 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 305 994 307 30 "" NOTE_INSN_DELETED_LABEL 288)

(insn 307 305 308 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 124 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 308 307 309 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 123 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 309 308 310 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 123 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 310 309 996 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 124 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 996 310 313 30 "" NOTE_INSN_DELETED_LABEL 363)

(note/s 313 996 315 30 "" NOTE_INSN_DELETED_LABEL 289)

(insn 315 313 316 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 126 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 316 315 317 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 125 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 317 316 318 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 125 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 318 317 892 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 126 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 892 318 320 30 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 372)) -1 (nil)
    (nil))
;; End of basic block 30, registers live:
 (nil)

(barrier 320 892 1008)

;; Start of basic block 33, registers live: (nil)
(code_label/s 1008 320 1011 33 366 "" [1 uses])

(note 1011 1008 1009 33 [bb 33] NOTE_INSN_BASIC_BLOCK)

(insn 1009 1011 1010 33 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1010 1009 321 33 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 321 1010 323 33 "" NOTE_INSN_DELETED_LABEL 290)

(insn 323 321 324 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 130 [ save_filt.192 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 324 323 325 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 129 [ save_eptr.191 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 325 324 326 33 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 326 325 327 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 115 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -32 [0xffffffffffffffe0])) [11 D.40971._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 327 326 328 33 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 328 327 329 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 115 [ this ]))) -1 (nil)
    (nil))

(jump_insn 329 328 331 33 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 358)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 33, registers live:
 (nil)

;; Start of basic block 35, registers live: (nil)
(note 331 329 332 35 [bb 35] NOTE_INSN_BASIC_BLOCK)

(note 332 331 333 35 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 333 332 334 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 214)
                (plus:DI (reg:DI 115 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 334 333 335 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 335 334 336 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 214)) -1 (nil)
    (nil))

(call_insn 336 335 842 35 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 36 [0x24])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 35, registers live:
 (nil)

;; Start of basic block 36, registers live: (nil)
(note 842 336 337 36 [bb 36] NOTE_INSN_BASIC_BLOCK)

(insn 337 842 339 36 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 116 [ D.41115 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 339 337 340 36 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 116 [ D.41115 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 340 339 342 36 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 358)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 36, registers live:
 (nil)

;; Start of basic block 37, registers live: (nil)
(note 342 340 343 37 [bb 37] NOTE_INSN_BASIC_BLOCK)

(note 343 342 344 37 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 344 343 345 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 215)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -6 [0xfffffffffffffffa])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 345 344 346 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 215)) -1 (nil)
    (nil))

(insn 346 345 347 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 115 [ this ])) -1 (nil)
    (nil))

(call_insn 347 346 348 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 348 347 349 37 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 358)) -1 (nil)
    (nil))
;; End of basic block 37, registers live:
 (nil)

(barrier 349 348 984)

;; Start of basic block 38, registers live: (nil)
(code_label/s 984 349 987 38 360 "" [1 uses])

(note 987 984 985 38 [bb 38] NOTE_INSN_BASIC_BLOCK)

(insn 985 987 986 38 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 986 985 350 38 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 350 986 352 38 "" NOTE_INSN_DELETED_LABEL 294)

(insn 352 350 353 38 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 118 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 353 352 354 38 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 117 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 354 353 355 38 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 117 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 355 354 886 38 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 118 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 886 355 357 38 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 364)) -1 (nil)
    (nil))
;; End of basic block 38, registers live:
 (nil)

(barrier 357 886 358)

;; Start of basic block 40, registers live: (nil)
(code_label 358 357 359 40 291 "" [3 uses])

(note 359 358 360 40 [bb 40] NOTE_INSN_BASIC_BLOCK)

(insn 360 359 361 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 129 [ save_eptr.191 ])) -1 (nil)
    (nil))

(insn 361 360 898 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 130 [ save_filt.192 ]))) -1 (nil)
    (nil))

(jump_insn 898 361 363 40 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 372)) -1 (nil)
    (nil))
;; End of basic block 40, registers live:
 (nil)

(barrier 363 898 988)

(note/s 988 363 364 129 "" NOTE_INSN_DELETED_LABEL 361)

;; Start of basic block 41, registers live: (nil)
(code_label/s 364 988 365 41 295 "" [2 uses])

(note 365 364 366 41 [bb 41] NOTE_INSN_BASIC_BLOCK)

(insn 366 365 367 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 120 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 367 366 368 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 119 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 368 367 369 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 119 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 369 368 888 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 120 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 888 369 371 41 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 372)) -1 (nil)
    (nil))
;; End of basic block 41, registers live:
 (nil)

(barrier 371 888 1012)

;; Start of basic block 42, registers live: (nil)
(code_label/s 1012 371 1015 42 367 "" [1 uses])

(note 1015 1012 1013 42 [bb 42] NOTE_INSN_BASIC_BLOCK)

(insn 1013 1015 1014 42 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1014 1013 372 42 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 42, registers live:
 (nil)

;; Start of basic block 43, registers live: (nil)
(code_label/s 372 1014 373 43 296 "" [4 uses])

(note 373 372 374 43 [bb 43] NOTE_INSN_BASIC_BLOCK)

(insn 374 373 375 43 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 132 [ save_filt.194 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 375 374 376 43 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 131 [ save_eptr.193 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 376 375 377 43 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 131 [ save_eptr.193 ])) -1 (nil)
    (nil))

(insn 377 376 900 43 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 132 [ save_filt.194 ]))) -1 (nil)
    (nil))

(jump_insn 900 377 379 43 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 510)) -1 (nil)
    (nil))
;; End of basic block 43, registers live:
 (nil)

(barrier 379 900 380)

;; Start of basic block 44, registers live: (nil)
(code_label 380 379 381 44 283 "" [1 uses])

(note 381 380 382 44 [bb 44] NOTE_INSN_BASIC_BLOCK)

(note 382 381 383 44 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 137)

(insn 383 382 384 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (parallel [
            (set (reg:DI 216)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -9 [0xfffffffffffffff7])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 384 383 385 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (set (reg:DI 217)
        (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -360 [0xfffffffffffffe98])) [11 s._M_stringbuf.D.34785._M_out_beg+0 S8 A64])) -1 (nil)
    (nil))

(insn 385 384 386 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (parallel [
            (set (reg:DI 218)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 386 385 387 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (set (reg:DI 2 cx)
        (reg:DI 216)) -1 (nil)
    (nil))

(insn 387 386 388 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (set (reg:DI 1 dx)
        (reg:DI 128 [ D.40961 ])) -1 (nil)
    (nil))

(insn 388 387 389 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (set (reg:DI 4 si)
        (reg:DI 217)) -1 (nil)
    (nil))

(insn 389 388 390 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (set (reg:DI 5 di)
        (reg:DI 218)) -1 (nil)
    (nil))

(call_insn 390 389 391 44 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:137 (call (mem:QI (symbol_ref/i:DI ("_ZNSsC1IPcEET_S1_RKSaIcE") [flags 0x1] <function_decl 0x2b4c43819d00 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 24 [0x18])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))
;; End of basic block 44, registers live:
 (nil)

;; Start of basic block 45, registers live: (nil)
(note 391 390 392 45 [bb 45] NOTE_INSN_BASIC_BLOCK)

(note 392 391 393 45 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 393 392 394 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 219)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 394 393 395 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 219)) -1 (nil)
    (nil))

(insn 395 394 396 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 98 [ this ])) -1 (nil)
    (nil))

(call_insn 396 395 843 45 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 25 [0x19])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 45, registers live:
 (nil)

;; Start of basic block 46, registers live: (nil)
(note 843 396 397 46 [bb 46] NOTE_INSN_BASIC_BLOCK)

(insn 397 843 399 46 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 140 [ D.40989 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 399 397 400 46 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 400 399 401 46 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 109 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -48 [0xffffffffffffffd0])) [11 D.40987._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 401 400 402 46 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 402 401 403 46 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 109 [ this ]))) -1 (nil)
    (nil))

(jump_insn 403 402 405 46 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 696)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 46, registers live:
 (nil)

;; Start of basic block 47, registers live: (nil)
(note 405 403 406 47 [bb 47] NOTE_INSN_BASIC_BLOCK)

(note 406 405 407 47 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 407 406 408 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 220)
                (plus:DI (reg:DI 109 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 408 407 409 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 409 408 410 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 220)) -1 (nil)
    (nil))

(call_insn 410 409 844 47 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 40 [0x28])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 47, registers live:
 (nil)

;; Start of basic block 48, registers live: (nil)
(note 844 410 411 48 [bb 48] NOTE_INSN_BASIC_BLOCK)

(insn 411 844 413 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 110 [ D.41171 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 413 411 414 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 110 [ D.41171 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 414 413 416 48 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 696)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 48, registers live:
 (nil)

;; Start of basic block 49, registers live: (nil)
(note 416 414 417 49 [bb 49] NOTE_INSN_BASIC_BLOCK)

(note 417 416 418 49 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 418 417 419 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 221)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -5 [0xfffffffffffffffb])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 419 418 420 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 221)) -1 (nil)
    (nil))

(insn 420 419 421 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 109 [ this ])) -1 (nil)
    (nil))

(call_insn 421 420 422 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 422 421 423 49 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 696)) -1 (nil)
    (nil))
;; End of basic block 49, registers live:
 (nil)

(barrier 423 422 976)

;; Start of basic block 50, registers live: (nil)
(code_label/s 976 423 979 50 358 "" [1 uses])

(note 979 976 977 50 [bb 50] NOTE_INSN_BASIC_BLOCK)

(insn 977 979 978 50 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 978 977 424 50 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 424 978 426 50 "" NOTE_INSN_DELETED_LABEL 299)

(insn 426 424 427 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 112 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 427 426 428 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 111 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 428 427 429 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 111 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 429 428 980 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 112 [ save_filt.160 ]))) -1 (nil)
    (nil))

(note/s 980 429 432 50 "" NOTE_INSN_DELETED_LABEL 359)

(note/s 432 980 434 50 "" NOTE_INSN_DELETED_LABEL 300)

(insn 434 432 435 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 114 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 435 434 436 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 113 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 436 435 437 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 113 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 437 436 884 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 114 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 884 437 439 50 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 491)) -1 (nil)
    (nil))
;; End of basic block 50, registers live:
 (nil)

(barrier 439 884 1000)

;; Start of basic block 53, registers live: (nil)
(code_label/s 1000 439 1003 53 364 "" [1 uses])

(note 1003 1000 1001 53 [bb 53] NOTE_INSN_BASIC_BLOCK)

(insn 1001 1003 1002 53 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1002 1001 440 53 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 440 1002 442 53 "" NOTE_INSN_DELETED_LABEL 301)

(insn 442 440 443 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 134 [ save_filt.196 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 443 442 444 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 133 [ save_eptr.195 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 444 443 445 53 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 445 444 446 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 103 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -48 [0xffffffffffffffd0])) [11 D.40987._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 446 445 447 53 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 447 446 448 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 103 [ this ]))) -1 (nil)
    (nil))

(jump_insn 448 447 450 53 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 477)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 53, registers live:
 (nil)

;; Start of basic block 55, registers live: (nil)
(note 450 448 451 55 [bb 55] NOTE_INSN_BASIC_BLOCK)

(note 451 450 452 55 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 452 451 453 55 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 222)
                (plus:DI (reg:DI 103 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 453 452 454 55 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 454 453 455 55 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 222)) -1 (nil)
    (nil))

(call_insn 455 454 845 55 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 44 [0x2c])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 55, registers live:
 (nil)

;; Start of basic block 56, registers live: (nil)
(note 845 455 456 56 [bb 56] NOTE_INSN_BASIC_BLOCK)

(insn 456 845 458 56 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 104 [ D.41227 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 458 456 459 56 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 104 [ D.41227 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 459 458 461 56 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 477)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 56, registers live:
 (nil)

;; Start of basic block 57, registers live: (nil)
(note 461 459 462 57 [bb 57] NOTE_INSN_BASIC_BLOCK)

(note 462 461 463 57 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 463 462 464 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 223)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -4 [0xfffffffffffffffc])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 464 463 465 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 223)) -1 (nil)
    (nil))

(insn 465 464 466 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 103 [ this ])) -1 (nil)
    (nil))

(call_insn 466 465 467 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 467 466 468 57 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 477)) -1 (nil)
    (nil))
;; End of basic block 57, registers live:
 (nil)

(barrier 468 467 968)

;; Start of basic block 58, registers live: (nil)
(code_label/s 968 468 971 58 356 "" [1 uses])

(note 971 968 969 58 [bb 58] NOTE_INSN_BASIC_BLOCK)

(insn 969 971 970 58 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 970 969 469 58 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 469 970 471 58 "" NOTE_INSN_DELETED_LABEL 305)

(insn 471 469 472 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 106 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 472 471 473 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 105 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 473 472 474 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 105 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 474 473 878 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 106 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 878 474 476 58 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 483)) -1 (nil)
    (nil))
;; End of basic block 58, registers live:
 (nil)

(barrier 476 878 477)

;; Start of basic block 60, registers live: (nil)
(code_label 477 476 478 60 302 "" [3 uses])

(note 478 477 479 60 [bb 60] NOTE_INSN_BASIC_BLOCK)

(insn 479 478 480 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 133 [ save_eptr.195 ])) -1 (nil)
    (nil))

(insn 480 479 894 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 134 [ save_filt.196 ]))) -1 (nil)
    (nil))

(jump_insn 894 480 482 60 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 491)) -1 (nil)
    (nil))
;; End of basic block 60, registers live:
 (nil)

(barrier 482 894 972)

(note/s 972 482 483 125 "" NOTE_INSN_DELETED_LABEL 357)

;; Start of basic block 61, registers live: (nil)
(code_label/s 483 972 484 61 306 "" [2 uses])

(note 484 483 485 61 [bb 61] NOTE_INSN_BASIC_BLOCK)

(insn 485 484 486 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 108 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 486 485 487 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 107 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 487 486 488 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 107 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 488 487 880 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 108 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 880 488 490 61 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 491)) -1 (nil)
    (nil))
;; End of basic block 61, registers live:
 (nil)

(barrier 490 880 1004)

;; Start of basic block 62, registers live: (nil)
(code_label/s 1004 490 1007 62 365 "" [1 uses])

(note 1007 1004 1005 62 [bb 62] NOTE_INSN_BASIC_BLOCK)

(insn 1005 1007 1006 62 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1006 1005 491 62 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 62, registers live:
 (nil)

;; Start of basic block 63, registers live: (nil)
(code_label/s 491 1006 492 63 307 "" [4 uses])

(note 492 491 493 63 [bb 63] NOTE_INSN_BASIC_BLOCK)

(insn 493 492 494 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 136 [ save_filt.198 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 494 493 495 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 135 [ save_eptr.197 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 495 494 496 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 135 [ save_eptr.197 ])) -1 (nil)
    (nil))

(insn 496 495 896 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 136 [ save_filt.198 ]))) -1 (nil)
    (nil))

(jump_insn 896 496 498 63 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 510)) -1 (nil)
    (nil))
;; End of basic block 63, registers live:
 (nil)

(barrier 498 896 499)

;; Start of basic block 64, registers live: (nil)
(code_label 499 498 500 64 281 "" [1 uses])

(note 500 499 501 64 [bb 64] NOTE_INSN_BASIC_BLOCK)

(note 501 500 502 64 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 486)

(insn 502 501 503 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 224)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 503 502 504 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (parallel [
            (set (reg:DI 225)
                (plus:DI (reg:DI 224)
                    (const_int 80 [0x50])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 504 503 505 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 4 si)
        (reg:DI 225)) -1 (nil)
    (nil))

(insn 505 504 506 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 5 di)
        (reg/f:DI 98 [ this ])) -1 (nil)
    (nil))

(call_insn 506 505 846 64 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSs6assignERKSs") [flags 0x41] <function_decl 0x2b4c42ae1000 assign>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 19 [0x13])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))
;; End of basic block 64, registers live:
 (nil)

;; Start of basic block 65, registers live: (nil)
(note 846 506 507 65 [bb 65] NOTE_INSN_BASIC_BLOCK)

(insn 507 846 508 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 141 [ D.40999 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(jump_insn 508 507 509 65 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (pc)
        (label_ref 696)) -1 (nil)
    (nil))
;; End of basic block 65, registers live:
 (nil)

(barrier 509 508 1016)

;; Start of basic block 66, registers live: (nil)
(code_label/s 1016 509 1019 66 368 "" [1 uses])

(note 1019 1016 1017 66 [bb 66] NOTE_INSN_BASIC_BLOCK)

(insn 1017 1019 1018 66 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1018 1017 510 66 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 66, registers live:
 (nil)

;; Start of basic block 67, registers live: (nil)
(code_label/s 510 1018 511 67 308 "" [3 uses])

(note 511 510 512 67 [bb 67] NOTE_INSN_BASIC_BLOCK)

(insn 512 511 513 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:SI 138 [ save_filt.200 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 513 512 514 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:486 (set (reg:DI 137 [ save_eptr.199 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 514 513 515 67 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 515 514 516 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 96 [ this ])
                (plus:DI (mem/s/f:DI (reg/f:DI 98 [ this ]) [11 <variable>._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 516 515 517 67 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 517 516 518 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 96 [ this ]))) -1 (nil)
    (nil))

(jump_insn 518 517 520 67 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 547)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 67, registers live:
 (nil)

;; Start of basic block 68, registers live: (nil)
(note 520 518 521 68 [bb 68] NOTE_INSN_BASIC_BLOCK)

(note 521 520 522 68 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 522 521 523 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 226)
                (plus:DI (reg:DI 96 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 523 522 524 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 524 523 525 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 226)) -1 (nil)
    (nil))

(call_insn 525 524 847 68 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 48 [0x30])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 68, registers live:
 (nil)

;; Start of basic block 69, registers live: (nil)
(note 847 525 526 69 [bb 69] NOTE_INSN_BASIC_BLOCK)

(insn 526 847 528 69 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 97 [ D.41284 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 528 526 529 69 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 97 [ D.41284 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 529 528 531 69 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 547)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 69, registers live:
 (nil)

;; Start of basic block 70, registers live: (nil)
(note 531 529 532 70 [bb 70] NOTE_INSN_BASIC_BLOCK)

(note 532 531 533 70 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 533 532 534 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 227)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -3 [0xfffffffffffffffd])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 534 533 535 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 227)) -1 (nil)
    (nil))

(insn 535 534 536 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 96 [ this ])) -1 (nil)
    (nil))

(call_insn 536 535 537 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 537 536 538 70 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 547)) -1 (nil)
    (nil))
;; End of basic block 70, registers live:
 (nil)

(barrier 538 537 960)

;; Start of basic block 71, registers live: (nil)
(code_label/s 960 538 963 71 354 "" [1 uses])

(note 963 960 961 71 [bb 71] NOTE_INSN_BASIC_BLOCK)

(insn 961 963 962 71 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 962 961 539 71 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 539 962 541 71 "" NOTE_INSN_DELETED_LABEL 312)

(insn 541 539 542 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 100 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 542 541 543 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 99 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 543 542 544 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 99 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 544 543 874 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 100 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 874 544 546 71 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 553)) -1 (nil)
    (nil))
;; End of basic block 71, registers live:
 (nil)

(barrier 546 874 547)

;; Start of basic block 73, registers live: (nil)
(code_label 547 546 548 73 309 "" [3 uses])

(note 548 547 549 73 [bb 73] NOTE_INSN_BASIC_BLOCK)

(insn 549 548 550 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 137 [ save_eptr.199 ])) -1 (nil)
    (nil))

(insn 550 549 902 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 138 [ save_filt.200 ]))) -1 (nil)
    (nil))

(jump_insn 902 550 552 73 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 561)) -1 (nil)
    (nil))
;; End of basic block 73, registers live:
 (nil)

(barrier 552 902 964)

(note/s 964 552 553 123 "" NOTE_INSN_DELETED_LABEL 355)

;; Start of basic block 74, registers live: (nil)
(code_label/s 553 964 554 74 313 "" [2 uses])

(note 554 553 555 74 [bb 74] NOTE_INSN_BASIC_BLOCK)

(insn 555 554 556 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 102 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 556 555 557 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 101 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 557 556 558 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 101 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 558 557 876 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 102 [ save_filt.162 ]))) -1 (nil)
    (nil))

(jump_insn 876 558 560 74 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 561)) -1 (nil)
    (nil))
;; End of basic block 74, registers live:
 (nil)

(barrier 560 876 1044)

;; Start of basic block 75, registers live: (nil)
(code_label/s 1044 560 1047 75 375 "" [1 uses])

(note 1047 1044 1045 75 [bb 75] NOTE_INSN_BASIC_BLOCK)

(insn 1045 1047 1046 75 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 1046 1045 561 75 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))
;; End of basic block 75, registers live:
 (nil)

;; Start of basic block 76, registers live: (nil)
(code_label/s 561 1046 562 76 314 "" [3 uses])

(note 562 561 563 76 [bb 76] NOTE_INSN_BASIC_BLOCK)

(insn 563 562 564 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 173 [ save_filt.178 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 564 563 565 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 174 [ save_eptr.177 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 565 564 566 76 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

(insn 566 565 567 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -400 [0xfffffffffffffe70])) [6 s.D.34935._vptr.basic_ostream+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436cad10 _ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 24 [0x18])))) -1 (nil)
    (nil))

(insn 567 566 568 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (mem/s:DI (reg:DI 160 [ D.40620 ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436cad10 _ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 64 [0x40])))) -1 (nil)
    (nil))

(note 568 567 569 76 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 569 568 570 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -392 [0xfffffffffffffe78])) [6 s._M_stringbuf.D.34785._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436ba370 _ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 570 569 571 76 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 571 570 572 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 84 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -320 [0xfffffffffffffec0])) [11 s._M_stringbuf._M_string._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 572 571 573 76 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 573 572 574 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 84 [ this ]))) -1 (nil)
    (nil))

(jump_insn 574 573 576 76 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 603)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 76, registers live:
 (nil)

;; Start of basic block 77, registers live: (nil)
(note 576 574 577 77 [bb 77] NOTE_INSN_BASIC_BLOCK)

(note 577 576 578 77 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 578 577 579 77 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 228)
                (plus:DI (reg:DI 84 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 579 578 580 77 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 580 579 581 77 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 228)) -1 (nil)
    (nil))

(call_insn 581 580 848 77 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 58 [0x3a])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 77, registers live:
 (nil)

;; Start of basic block 78, registers live: (nil)
(note 848 581 582 78 [bb 78] NOTE_INSN_BASIC_BLOCK)

(insn 582 848 584 78 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 85 [ D.41393 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 584 582 585 78 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 85 [ D.41393 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 585 584 587 78 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 603)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 78, registers live:
 (nil)

;; Start of basic block 79, registers live: (nil)
(note 587 585 588 79 [bb 79] NOTE_INSN_BASIC_BLOCK)

(note 588 587 589 79 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 589 588 590 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 229)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -2 [0xfffffffffffffffe])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 590 589 591 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 229)) -1 (nil)
    (nil))

(insn 591 590 592 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 84 [ this ])) -1 (nil)
    (nil))

(call_insn 592 591 593 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 593 592 594 79 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 603)) -1 (nil)
    (nil))
;; End of basic block 79, registers live:
 (nil)

(barrier 594 593 940)

;; Start of basic block 80, registers live: (nil)
(code_label/s 940 594 943 80 349 "" [1 uses])

(note 943 940 941 80 [bb 80] NOTE_INSN_BASIC_BLOCK)

(insn 941 943 942 80 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 942 941 595 80 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 595 942 597 80 "" NOTE_INSN_DELETED_LABEL 318)

(insn 597 595 598 80 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 87 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 598 597 599 80 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 86 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 599 598 600 80 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 86 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 600 599 863 80 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 87 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 863 600 602 80 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 615)) -1 (nil)
    (nil))
;; End of basic block 80, registers live:
 (nil)

(barrier 602 863 603)

;; Start of basic block 82, registers live: (nil)
(code_label 603 602 604 82 315 "" [3 uses])

(note 604 603 605 82 [bb 82] NOTE_INSN_BASIC_BLOCK)

(note 605 604 606 82 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 606 605 607 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 230)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 607 606 608 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 82 [ this ])
                (plus:DI (reg:DI 230)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 608 607 609 82 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 200)

(insn 609 608 610 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (mem/s:DI (reg:DI 82 [ this ]) [6 <variable>._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 610 609 611 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (parallel [
            (set (reg:DI 231)
                (plus:DI (reg:DI 82 [ this ])
                    (const_int 56 [0x38])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 611 610 612 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 5 di)
        (reg:DI 231)) -1 (nil)
    (nil))

(call_insn 612 611 613 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeD1Ev") [flags 0x41] <function_decl 0x2b4c42c1d300 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 613 612 614 82 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (pc)
        (label_ref 642)) -1 (nil)
    (nil))
;; End of basic block 82, registers live:
 (nil)

(barrier 614 613 944)

(note/s 944 614 615 118 "" NOTE_INSN_DELETED_LABEL 350)

;; Start of basic block 83, registers live: (nil)
(code_label/s 615 944 616 83 320 "" [2 uses])

(note 616 615 617 83 [bb 83] NOTE_INSN_BASIC_BLOCK)

(insn 617 616 618 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:SI 89 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 618 617 619 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 88 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 619 618 620 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 196)
        (reg:DI 88 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 620 619 948 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 89 [ save_filt.162 ]))) -1 (nil)
    (nil))

(note/s 948 620 623 83 "" NOTE_INSN_DELETED_LABEL 351)

(note/s 623 948 625 83 "" NOTE_INSN_DELETED_LABEL 321)

(insn 625 623 626 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:SI 91 [ save_filt.190 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 626 625 627 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 90 [ save_eptr.189 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 627 626 628 83 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 628 627 629 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 232)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 629 628 630 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 83 [ this ])
                (plus:DI (reg:DI 232)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 630 629 631 83 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 200)

(insn 631 630 632 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (mem/s:DI (reg:DI 83 [ this ]) [6 <variable>._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 632 631 633 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (parallel [
            (set (reg:DI 233)
                (plus:DI (reg:DI 83 [ this ])
                    (const_int 56 [0x38])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 633 632 634 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 5 di)
        (reg:DI 233)) -1 (nil)
    (nil))

(call_insn 634 633 637 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeD1Ev") [flags 0x41] <function_decl 0x2b4c42c1d300 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 637 634 638 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 196)
        (reg:DI 90 [ save_eptr.189 ])) -1 (nil)
    (nil))

(insn 638 637 867 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 91 [ save_filt.190 ]))) -1 (nil)
    (nil))

(jump_insn 867 638 640 83 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (pc)
        (label_ref 652)) -1 (nil)
    (nil))
;; End of basic block 83, registers live:
 (nil)

(barrier 640 867 641)

(note 641 640 642 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

;; Start of basic block 85, registers live: (nil)
(code_label 642 641 643 85 319 "" [1 uses])

(note 643 642 644 85 [bb 85] NOTE_INSN_BASIC_BLOCK)

(note 644 643 645 85 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

(insn 645 644 646 85 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (parallel [
            (set (reg:DI 80 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 646 645 647 85 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream") 113)

(insn 647 646 648 85 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (reg:DI 80 [ this ]) [6 <variable>._vptr.basic_ostream+0 S8 A64])
        (reg:DI 169 [ D.40609 ])) -1 (nil)
    (nil))

(insn 648 647 649 85 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 234)
        (mem:DI (reg:DI 168 [ D.40611 ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 649 648 650 85 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (plus:DI (reg:DI 80 [ this ])
                (reg:DI 234)) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (reg:DI 167 [ iftmp.131 ])) -1 (nil)
    (nil))

(jump_insn 650 649 651 85 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (pc)
        (label_ref 668)) -1 (nil)
    (nil))
;; End of basic block 85, registers live:
 (nil)

(barrier 651 650 952)

(note/s 952 651 652 120 "" NOTE_INSN_DELETED_LABEL 352)

;; Start of basic block 86, registers live: (nil)
(code_label/s 652 952 653 86 324 "" [2 uses])

(note 653 652 654 86 [bb 86] NOTE_INSN_BASIC_BLOCK)

(insn 654 653 655 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:SI 93 [ save_filt.204 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 655 654 656 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 92 [ save_eptr.203 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 656 655 657 86 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

(insn 657 656 658 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (parallel [
            (set (reg:DI 81 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 658 657 659 86 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream") 113)

(insn 659 658 660 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (reg:DI 81 [ this ]) [6 <variable>._vptr.basic_ostream+0 S8 A64])
        (reg:DI 169 [ D.40609 ])) -1 (nil)
    (nil))

(insn 660 659 661 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 235)
        (mem:DI (reg:DI 168 [ D.40611 ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 661 660 664 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (plus:DI (reg:DI 81 [ this ])
                (reg:DI 235)) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (reg:DI 167 [ iftmp.131 ])) -1 (nil)
    (nil))

(insn 664 661 665 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 196)
        (reg:DI 92 [ save_eptr.203 ])) -1 (nil)
    (nil))

(insn 665 664 869 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 93 [ save_filt.204 ]))) -1 (nil)
    (nil))

(jump_insn 869 665 667 86 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (pc)
        (label_ref 676)) -1 (nil)
    (nil))
;; End of basic block 86, registers live:
 (nil)

(barrier 667 869 668)

;; Start of basic block 87, registers live: (nil)
(code_label 668 667 669 87 323 "" [1 uses])

(note 669 668 670 87 [bb 87] NOTE_INSN_BASIC_BLOCK)

(note 670 669 671 87 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h") 268)

(insn 671 670 672 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (mem/s:DI (reg:DI 157 [ this ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9basic_iosIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c430446e0 _ZTVSt9basic_iosIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 672 671 673 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 157 [ this ])) -1 (nil)
    (nil))

(call_insn 673 672 674 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_baseD2Ev") [flags 0x41] <function_decl 0x2b4c42c93900 __base_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 2 [0x2])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 674 673 675 87 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (pc)
        (label_ref 690)) -1 (nil)
    (nil))
;; End of basic block 87, registers live:
 (nil)

(barrier 675 674 956)

(note/s 956 675 676 121 "" NOTE_INSN_DELETED_LABEL 353)

;; Start of basic block 88, registers live: (nil)
(code_label/s 676 956 677 88 327 "" [2 uses])

(note 677 676 678 88 [bb 88] NOTE_INSN_BASIC_BLOCK)

(insn 678 677 679 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:SI 95 [ save_filt.206 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 679 678 680 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 94 [ save_eptr.205 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 680 679 681 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (mem/s:DI (reg:DI 157 [ this ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9basic_iosIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c430446e0 _ZTVSt9basic_iosIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 681 680 682 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 157 [ this ])) -1 (nil)
    (nil))

(call_insn 682 681 685 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_baseD2Ev") [flags 0x41] <function_decl 0x2b4c42c93900 __base_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 54 [0x36])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 685 682 686 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 196)
        (reg:DI 94 [ save_eptr.205 ])) -1 (nil)
    (nil))

(insn 686 685 871 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 95 [ save_filt.206 ]))) -1 (nil)
    (nil))

(insn 871 686 872 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 196)) -1 (nil)
    (nil))

(call_insn 872 871 688 88 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 88, registers live:
 (nil)

(barrier 688 872 689)

(note 689 688 690 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

;; Start of basic block 89, registers live: (nil)
(code_label 690 689 691 89 326 "" [1 uses])

(note 691 690 692 89 [bb 89] NOTE_INSN_BASIC_BLOCK)

(insn 692 691 693 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (reg:DI 196)
        (reg:DI 174 [ save_eptr.177 ])) -1 (nil)
    (nil))

(insn 693 692 917 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 173 [ save_filt.178 ]))) -1 (nil)
    (nil))

(insn 917 693 918 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (reg:DI 5 di)
        (reg:DI 196)) -1 (nil)
    (nil))

(call_insn 918 917 695 89 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 89, registers live:
 (nil)

(barrier 695 918 696)

;; Start of basic block 90, registers live: (nil)
(code_label 696 695 697 90 285 "" [7 uses])

(note 697 696 698 90 [bb 90] NOTE_INSN_BASIC_BLOCK)

(insn 698 697 699 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -400 [0xfffffffffffffe70])) [6 s.D.34935._vptr.basic_ostream+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436cad10 _ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 24 [0x18])))) -1 (nil)
    (nil))

(insn 699 698 700 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (set (mem/s:DI (reg:DI 160 [ D.40620 ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436cad10 _ZTVSt19basic_ostringstreamIcSt11char_traitsIcESaIcEE>)
                (const_int 64 [0x40])))) -1 (nil)
    (nil))

(note 700 699 701 90 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 701 700 702 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (set (mem/s/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -392 [0xfffffffffffffe78])) [6 s._M_stringbuf.D.34785._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE") [flags 0x40] <var_decl 0x2b4c436ba370 _ZTVSt15basic_stringbufIcSt11char_traitsIcESaIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(note 702 701 703 90 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 283)

(insn 703 702 704 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:283 (parallel [
            (set (reg:DI 68 [ this ])
                (plus:DI (mem/s/f/c:DI (plus:DI (reg/f:DI 54 virtual-stack-vars)
                            (const_int -320 [0xfffffffffffffec0])) [11 s._M_stringbuf._M_string._M_dataplus._M_p+0 S8 A64])
                    (const_int -24 [0xffffffffffffffe8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 704 703 705 90 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 230)

(insn 705 704 706 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 156 [ this ])
            (reg:DI 68 [ this ]))) -1 (nil)
    (nil))

(jump_insn 706 705 708 90 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:230 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 735)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9901 [0x26ad])
        (expr_list:REG_BR_PRED (concat (const_int 6 [0x6])
                (const_int 9901 [0x26ad]))
            (nil))))
;; End of basic block 90, registers live:
 (nil)

;; Start of basic block 91, registers live: (nil)
(note 708 706 709 91 [bb 91] NOTE_INSN_BASIC_BLOCK)

(note 709 708 710 91 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 232)

(insn 710 709 711 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (parallel [
            (set (reg:DI 236)
                (plus:DI (reg:DI 68 [ this ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 711 710 712 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 4 si)
        (const_int -1 [0xffffffffffffffff])) -1 (nil)
    (nil))

(insn 712 711 713 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:DI 5 di)
        (reg:DI 236)) -1 (nil)
    (nil))

(call_insn 713 712 849 91 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZN9__gnu_cxx18__exchange_and_addEPVii") [flags 0x41] <function_decl 0x2b4c42588b00 __exchange_and_add>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 68 [0x44])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 91, registers live:
 (nil)

;; Start of basic block 92, registers live: (nil)
(note 849 713 714 92 [bb 92] NOTE_INSN_BASIC_BLOCK)

(insn 714 849 716 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:SI 69 [ D.41598 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 716 714 717 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg:SI 69 [ D.41598 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 717 716 719 92 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:232 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 735)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 8977 [0x2311])
        (nil)))
;; End of basic block 92, registers live:
 (nil)

;; Start of basic block 93, registers live: (nil)
(note 719 717 720 93 [bb 93] NOTE_INSN_BASIC_BLOCK)

(note 720 719 721 93 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h") 233)

(insn 721 720 722 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (parallel [
            (set (reg:DI 237)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 722 721 723 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 4 si)
        (reg:DI 237)) -1 (nil)
    (nil))

(insn 723 722 724 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 5 di)
        (reg:DI 68 [ this ])) -1 (nil)
    (nil))

(call_insn 724 723 725 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (call (mem:QI (symbol_ref:DI ("_ZNSs4_Rep10_M_destroyERKSaIcE") [flags 0x41] <function_decl 0x2b4c42b2eb00 _M_destroy>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(jump_insn 725 724 726 93 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 735)) -1 (nil)
    (nil))
;; End of basic block 93, registers live:
 (nil)

(barrier 726 725 920)

;; Start of basic block 94, registers live: (nil)
(code_label/s 920 726 923 94 344 "" [1 uses])

(note 923 920 921 94 [bb 94] NOTE_INSN_BASIC_BLOCK)

(insn 921 923 922 94 (set (reg:DI 196)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 922 921 727 94 (set (reg:DI 195)
        (reg:DI 1 dx)) -1 (nil)
    (nil))

(note/s 727 922 729 94 "" NOTE_INSN_DELETED_LABEL 332)

(insn 729 727 730 94 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:SI 71 [ save_filt.160 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 730 729 731 94 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 70 [ save_eptr.159 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 731 730 732 94 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 196)
        (reg:DI 70 [ save_eptr.159 ])) -1 (nil)
    (nil))

(insn 732 731 852 94 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 71 [ save_filt.160 ]))) -1 (nil)
    (nil))

(jump_insn 852 732 734 94 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_string.h:233 (set (pc)
        (label_ref 747)) -1 (nil)
    (nil))
;; End of basic block 94, registers live:
 (nil)

(barrier 734 852 735)

;; Start of basic block 96, registers live: (nil)
(code_label 735 734 736 96 329 "" [3 uses])

(note 736 735 737 96 [bb 96] NOTE_INSN_BASIC_BLOCK)

(note 737 736 738 96 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 738 737 739 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 238)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 739 738 740 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 66 [ this ])
                (plus:DI (reg:DI 238)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 740 739 741 96 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 200)

(insn 741 740 742 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (mem/s:DI (reg:DI 66 [ this ]) [6 <variable>._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 742 741 743 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (parallel [
            (set (reg:DI 239)
                (plus:DI (reg:DI 66 [ this ])
                    (const_int 56 [0x38])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 743 742 744 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 5 di)
        (reg:DI 239)) -1 (nil)
    (nil))

(call_insn 744 743 745 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeD1Ev") [flags 0x41] <function_decl 0x2b4c42c1d300 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 745 744 746 96 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (pc)
        (label_ref 774)) -1 (nil)
    (nil))
;; End of basic block 96, registers live:
 (nil)

(barrier 746 745 924)

(note/s 924 746 747 113 "" NOTE_INSN_DELETED_LABEL 345)

;; Start of basic block 97, registers live: (nil)
(code_label/s 747 924 748 97 334 "" [2 uses])

(note 748 747 749 97 [bb 97] NOTE_INSN_BASIC_BLOCK)

(insn 749 748 750 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:SI 73 [ save_filt.162 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 750 749 751 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 72 [ save_eptr.161 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 751 750 752 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 196)
        (reg:DI 72 [ save_eptr.161 ])) -1 (nil)
    (nil))

(insn 752 751 928 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 73 [ save_filt.162 ]))) -1 (nil)
    (nil))

(note/s 928 752 755 97 "" NOTE_INSN_DELETED_LABEL 346)

(note/s 755 928 757 97 "" NOTE_INSN_DELETED_LABEL 335)

(insn 757 755 758 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:SI 75 [ save_filt.190 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 758 757 759 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 74 [ save_eptr.189 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 759 758 760 97 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

(insn 760 759 761 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 240)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 761 760 762 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd:71 (parallel [
            (set (reg:DI 67 [ this ])
                (plus:DI (reg:DI 240)
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 762 761 763 97 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf") 200)

(insn 763 762 764 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (mem/s:DI (reg:DI 67 [ this ]) [6 <variable>._vptr.basic_streambuf+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt15basic_streambufIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c42d1f420 _ZTVSt15basic_streambufIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 764 763 765 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (parallel [
            (set (reg:DI 241)
                (plus:DI (reg:DI 67 [ this ])
                    (const_int 56 [0x38])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 765 764 766 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 5 di)
        (reg:DI 241)) -1 (nil)
    (nil))

(call_insn 766 765 769 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (call (mem:QI (symbol_ref:DI ("_ZNSt6localeD1Ev") [flags 0x41] <function_decl 0x2b4c42c1d300 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 769 766 770 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 196)
        (reg:DI 74 [ save_eptr.189 ])) -1 (nil)
    (nil))

(insn 770 769 856 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 75 [ save_filt.190 ]))) -1 (nil)
    (nil))

(jump_insn 856 770 772 97 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/streambuf:200 (set (pc)
        (label_ref 784)) -1 (nil)
    (nil))
;; End of basic block 97, registers live:
 (nil)

(barrier 772 856 773)

(note 773 772 774 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iosfwd") 71)

;; Start of basic block 99, registers live: (nil)
(code_label 774 773 775 99 333 "" [1 uses])

(note 775 774 776 99 [bb 99] NOTE_INSN_BASIC_BLOCK)

(note 776 775 777 99 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

(insn 777 776 778 99 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (parallel [
            (set (reg:DI 64 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 778 777 779 99 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream") 113)

(insn 779 778 780 99 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (reg:DI 64 [ this ]) [6 <variable>._vptr.basic_ostream+0 S8 A64])
        (reg:DI 169 [ D.40609 ])) -1 (nil)
    (nil))

(insn 780 779 781 99 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 242)
        (mem:DI (reg:DI 168 [ D.40611 ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 781 780 782 99 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (plus:DI (reg:DI 64 [ this ])
                (reg:DI 242)) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (reg:DI 167 [ iftmp.131 ])) -1 (nil)
    (nil))

(jump_insn 782 781 783 99 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (pc)
        (label_ref 800)) -1 (nil)
    (nil))
;; End of basic block 99, registers live:
 (nil)

(barrier 783 782 932)

(note/s 932 783 784 115 "" NOTE_INSN_DELETED_LABEL 347)

;; Start of basic block 100, registers live: (nil)
(code_label/s 784 932 785 100 338 "" [2 uses])

(note 785 784 786 100 [bb 100] NOTE_INSN_BASIC_BLOCK)

(insn 786 785 787 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:SI 77 [ save_filt.204 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 787 786 788 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 76 [ save_eptr.203 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(note 788 787 789 100 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

(insn 789 788 790 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream:444 (parallel [
            (set (reg:DI 65 [ this ])
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -400 [0xfffffffffffffe70])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 790 789 791 100 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream") 113)

(insn 791 790 792 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (reg:DI 65 [ this ]) [6 <variable>._vptr.basic_ostream+0 S8 A64])
        (reg:DI 169 [ D.40609 ])) -1 (nil)
    (nil))

(insn 792 791 793 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 243)
        (mem:DI (reg:DI 168 [ D.40611 ]) [19 S8 A64])) -1 (nil)
    (nil))

(insn 793 792 796 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (mem/s:DI (plus:DI (reg:DI 65 [ this ])
                (reg:DI 243)) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (reg:DI 167 [ iftmp.131 ])) -1 (nil)
    (nil))

(insn 796 793 797 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 196)
        (reg:DI 76 [ save_eptr.203 ])) -1 (nil)
    (nil))

(insn 797 796 858 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 77 [ save_filt.204 ]))) -1 (nil)
    (nil))

(jump_insn 858 797 799 100 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/ostream:113 (set (pc)
        (label_ref 808)) -1 (nil)
    (nil))
;; End of basic block 100, registers live:
 (nil)

(barrier 799 858 800)

;; Start of basic block 101, registers live: (nil)
(code_label 800 799 801 101 337 "" [1 uses])

(note 801 800 802 101 [bb 101] NOTE_INSN_BASIC_BLOCK)

(note 802 801 803 101 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h") 268)

(insn 803 802 804 101 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (mem/s:DI (reg:DI 157 [ this ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9basic_iosIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c430446e0 _ZTVSt9basic_iosIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 804 803 805 101 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 157 [ this ])) -1 (nil)
    (nil))

(call_insn 805 804 806 101 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_baseD2Ev") [flags 0x41] <function_decl 0x2b4c42c93900 __base_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 806 805 807 101 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (pc)
        (label_ref:DI 828)) 548 {jump} (nil)
    (nil))
;; End of basic block 101, registers live:
 (nil)

(barrier 807 806 936)

(note/s 936 807 808 116 "" NOTE_INSN_DELETED_LABEL 348)

;; Start of basic block 102, registers live: (nil)
(code_label/s 808 936 809 102 341 "" [2 uses])

(note 809 808 810 102 [bb 102] NOTE_INSN_BASIC_BLOCK)

(insn 810 809 811 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:SI 79 [ save_filt.206 ])
        (subreg:SI (reg:DI 195) 0)) -1 (nil)
    (nil))

(insn 811 810 812 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 78 [ save_eptr.205 ])
        (reg:DI 196)) -1 (nil)
    (nil))

(insn 812 811 813 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (mem/s:DI (reg:DI 157 [ this ]) [6 <variable>.D.25875._vptr.ios_base+0 S8 A64])
        (const:DI (plus:DI (symbol_ref:DI ("_ZTVSt9basic_iosIcSt11char_traitsIcEE") [flags 0x40] <var_decl 0x2b4c430446e0 _ZTVSt9basic_iosIcSt11char_traitsIcEE>)
                (const_int 16 [0x10])))) -1 (nil)
    (nil))

(insn 813 812 814 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 157 [ this ])) -1 (nil)
    (nil))

(call_insn 814 813 817 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_baseD2Ev") [flags 0x41] <function_decl 0x2b4c42c93900 __base_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 64 [0x40])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 817 814 818 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 196)
        (reg:DI 78 [ save_eptr.205 ])) -1 (nil)
    (nil))

(insn 818 817 860 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 195)
        (sign_extend:DI (reg:SI 79 [ save_filt.206 ]))) -1 (nil)
    (nil))

(insn 860 818 861 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (set (reg:DI 5 di)
        (reg:DI 196)) -1 (nil)
    (nil))

(call_insn 861 860 820 102 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/basic_ios.h:268 (call (mem:QI (symbol_ref:DI ("_Unwind_Resume") [flags 0x41]) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 102, registers live:
 (nil)

(barrier 820 861 821)

(note 821 820 826 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/sstream") 444)

(note 826 821 827 NOTE_INSN_FUNCTION_END)

(note 827 826 828 ("./CppUtilities.cc") 51)

;; Start of basic block 103, registers live: (nil)
(code_label 828 827 850 103 261 "" [1 uses])

(note 850 828 829 103 [bb 103] NOTE_INSN_BASIC_BLOCK)

(insn 829 850 830 103 ./CppUtilities.cc:51 (set (reg/i:DI 0 ax)
        (reg/f:DI 180 [ <result> ])) -1 (nil)
    (nil))

(insn 830 829 836 103 ./CppUtilities.cc:51 (set (reg/i:DI 0 ax)
        (reg/f:DI 180 [ <result> ])) -1 (nil)
    (nil))

(insn 836 830 0 103 ./CppUtilities.cc:51 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 103, registers live:
 (nil)

