
;; Function summaryStats_t::summaryStats_t() (_ZN14summaryStats_tC2Ev)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 3 failed.
Merged 0 and 2 without moving.
Merged 0 and 3 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 9 ("./CppStatUtilities.cc") 32)

;; Start of basic block 0, registers live: (nil)
(note 9 1 6 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 6 9 7 0 ./CppStatUtilities.cc:32 (set (reg/f:DI 58 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("./CppStatUtilities.cc") 33)

(insn 12 11 13 0 ./CppStatUtilities.cc:33 (set (reg:DF 59)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (reg/f:DI 58 [ this ]) [3 <variable>.min+0 S8 A64])
        (reg:DF 59)) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppStatUtilities.cc:33 (set (reg:DF 60)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])
        (reg:DF 60)) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppStatUtilities.cc:33 (set (reg:DF 61)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])
        (reg:DF 61)) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:33 (set (reg:DF 62)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 19 18 20 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])
        (reg:DF 62)) -1 (nil)
    (nil))

(insn 20 19 21 0 ./CppStatUtilities.cc:33 (set (reg:DF 63)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])
        (reg:DF 63)) -1 (nil)
    (nil))

(insn 22 21 23 0 ./CppStatUtilities.cc:33 (set (reg:DF 64)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 23 22 24 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])
        (reg:DF 64)) -1 (nil)
    (nil))

(insn 24 23 25 0 ./CppStatUtilities.cc:33 (set (reg:DF 65)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 25 24 26 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])
        (reg:DF 65)) -1 (nil)
    (nil))
;; End of basic block 0, registers live:
 (nil)

(note 26 25 27 NOTE_INSN_FUNCTION_END)

(note 27 26 0 ("./CppStatUtilities.cc") 34)


;; Function summaryStats_t::summaryStats_t() (_ZN14summaryStats_tC1Ev)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Forwarding edge 0->2 to 3 failed.
Merged 0 and 2 without moving.
Merged 0 and 3 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./CppStatUtilities.cc") 32)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./CppStatUtilities.cc:32 (set (reg/f:DI 58 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./CppStatUtilities.cc") 33)

(insn 9 8 10 0 ./CppStatUtilities.cc:33 (set (reg:DF 59)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (reg/f:DI 58 [ this ]) [3 <variable>.min+0 S8 A64])
        (reg:DF 59)) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppStatUtilities.cc:33 (set (reg:DF 60)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])
        (reg:DF 60)) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:33 (set (reg:DF 61)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])
        (reg:DF 61)) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:33 (set (reg:DF 62)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])
        (reg:DF 62)) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:33 (set (reg:DF 63)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])
        (reg:DF 63)) -1 (nil)
    (nil))

(insn 19 18 20 0 ./CppStatUtilities.cc:33 (set (reg:DF 64)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 20 19 21 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])
        (reg:DF 64)) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppStatUtilities.cc:33 (set (reg:DF 65)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 22 21 23 0 ./CppStatUtilities.cc:33 (set (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])
        (reg:DF 65)) -1 (nil)
    (nil))
;; End of basic block 0, registers live:
 (nil)

(note 23 22 24 NOTE_INSN_FUNCTION_END)

(note 24 23 0 ("./CppStatUtilities.cc") 34)


;; Function summaryStats_t& summaryStats_t::operator=(const summaryStats_t&) (_ZN14summaryStats_taSERKS_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 3 and 4 without moving.
Merged 3 and 5 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("./CppStatUtilities.cc") 37)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppStatUtilities.cc:37 (set (reg/f:DI 59 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:37 (set (reg/v/f:DI 60 [ other ])
        (reg:DI 4 si [ other ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppStatUtilities.cc") 39)

(insn 10 9 11 0 ./CppStatUtilities.cc:39 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/f:DI 59 [ this ])
            (reg/v/f:DI 60 [ other ]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./CppStatUtilities.cc:39 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 35)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1036 [0x40c])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 13 11 14 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 14 13 15 1 ("./CppStatUtilities.cc") 41)

(insn 15 14 16 1 ./CppStatUtilities.cc:41 (set (reg:DF 61)
        (mem/s:DF (reg/v/f:DI 60 [ other ]) [3 <variable>.min+0 S8 A64])) -1 (nil)
    (nil))

(insn 16 15 17 1 ./CppStatUtilities.cc:41 (set (mem/s:DF (reg/f:DI 59 [ this ]) [3 <variable>.min+0 S8 A64])
        (reg:DF 61)) -1 (nil)
    (nil))

(note 17 16 18 1 ("./CppStatUtilities.cc") 42)

(insn 18 17 19 1 ./CppStatUtilities.cc:42 (set (reg:DF 62)
        (mem/s:DF (plus:DI (reg/v/f:DI 60 [ other ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])) -1 (nil)
    (nil))

(insn 19 18 20 1 ./CppStatUtilities.cc:42 (set (mem/s:DF (plus:DI (reg/f:DI 59 [ this ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])
        (reg:DF 62)) -1 (nil)
    (nil))

(note 20 19 21 1 ("./CppStatUtilities.cc") 43)

(insn 21 20 22 1 ./CppStatUtilities.cc:43 (set (reg:DF 63)
        (mem/s:DF (plus:DI (reg/v/f:DI 60 [ other ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])) -1 (nil)
    (nil))

(insn 22 21 23 1 ./CppStatUtilities.cc:43 (set (mem/s:DF (plus:DI (reg/f:DI 59 [ this ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])
        (reg:DF 63)) -1 (nil)
    (nil))

(note 23 22 24 1 ("./CppStatUtilities.cc") 44)

(insn 24 23 25 1 ./CppStatUtilities.cc:44 (set (reg:DF 64)
        (mem/s:DF (plus:DI (reg/v/f:DI 60 [ other ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])) -1 (nil)
    (nil))

(insn 25 24 26 1 ./CppStatUtilities.cc:44 (set (mem/s:DF (plus:DI (reg/f:DI 59 [ this ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])
        (reg:DF 64)) -1 (nil)
    (nil))

(note 26 25 27 1 ("./CppStatUtilities.cc") 45)

(insn 27 26 28 1 ./CppStatUtilities.cc:45 (set (reg:DF 65)
        (mem/s:DF (plus:DI (reg/v/f:DI 60 [ other ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])) -1 (nil)
    (nil))

(insn 28 27 29 1 ./CppStatUtilities.cc:45 (set (mem/s:DF (plus:DI (reg/f:DI 59 [ this ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])
        (reg:DF 65)) -1 (nil)
    (nil))

(note 29 28 30 1 ("./CppStatUtilities.cc") 46)

(insn 30 29 31 1 ./CppStatUtilities.cc:46 (set (reg:DF 66)
        (mem/s:DF (plus:DI (reg/v/f:DI 60 [ other ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])) -1 (nil)
    (nil))

(insn 31 30 32 1 ./CppStatUtilities.cc:46 (set (mem/s:DF (plus:DI (reg/f:DI 59 [ this ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])
        (reg:DF 66)) -1 (nil)
    (nil))

(note 32 31 33 1 ("./CppStatUtilities.cc") 47)

(insn 33 32 34 1 ./CppStatUtilities.cc:47 (set (reg:DF 67)
        (mem/s:DF (plus:DI (reg/v/f:DI 60 [ other ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])) -1 (nil)
    (nil))

(insn 34 33 35 1 ./CppStatUtilities.cc:47 (set (mem/s:DF (plus:DI (reg/f:DI 59 [ this ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])
        (reg:DF 67)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 35 34 36 2 6 "" [1 uses])

(note 36 35 37 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 37 36 40 2 ./CppStatUtilities.cc:47 (set (reg:DI 58 [ <result> ])
        (reg/f:DI 59 [ this ])) -1 (nil)
    (nil))

(note 40 37 41 2 NOTE_INSN_FUNCTION_END)

(note 41 40 43 2 ("./CppStatUtilities.cc") 51)

(insn 43 41 49 2 ./CppStatUtilities.cc:51 (set (reg/i:DI 0 ax)
        (reg:DI 58 [ <result> ])) -1 (nil)
    (nil))

(insn 49 43 0 2 ./CppStatUtilities.cc:51 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)


;; Function void whichMax(double*, int, double&, int&) (_Z8whichMaxPdiRdRi)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->8 redirected to 10
Deleted label in block 2.
Edge 3->5 redirected to 6
Forwarding edge 3->4 to 7 failed.
Forwarding edge 3->4 to 7 failed.
Deleting block 5.
Forwarding edge 7->8 to 10 failed.
Deleted label in block 8.
Deleting fallthru block 8.
Forwarding edge 7->9 to 10 failed.
Deleting fallthru block 9.


try_optimize_cfg iteration 2

Forwarding edge 3->4 to 7 failed.


try_optimize_cfg iteration 1

Forwarding edge 2->3 to 5 failed.
(note 1 0 11 ("./CppStatUtilities.cc") 134)

;; Start of basic block 0, registers live: (nil)
(note 11 1 5 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 5 11 6 0 ./CppStatUtilities.cc:134 (set (reg/v/f:DI 63 [ data ])
        (reg:DI 5 di [ data ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:134 (set (reg/v:SI 64 [ dataLen ])
        (reg:SI 4 si [ dataLen ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:134 (set (reg/v/f:DI 65 [ maxVal ])
        (reg:DI 1 dx [ maxVal ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppStatUtilities.cc:134 (set (reg/v/f:DI 66 [ maxInd ])
        (reg:DI 2 cx [ maxInd ])) -1 (nil)
    (nil))

(note 9 8 13 0 NOTE_INSN_FUNCTION_BEG)

(note 13 9 14 0 ("./CppStatUtilities.cc") 136)

(insn 14 13 15 0 ./CppStatUtilities.cc:136 (set (reg:DF 67)
        (mem:DF (reg/v/f:DI 63 [ data ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:136 (set (mem:DF (reg/v/f:DI 65 [ maxVal ]) [3 S8 A64])
        (reg:DF 67)) -1 (nil)
    (nil))

(note 16 15 17 0 ("./CppStatUtilities.cc") 137)

(insn 17 16 18 0 ./CppStatUtilities.cc:137 (set (mem:SI (reg/v/f:DI 66 [ maxInd ]) [5 S4 A32])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 18 17 19 0 ("./CppStatUtilities.cc") 140)

(insn 19 18 20 0 ./CppStatUtilities.cc:140 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 64 [ dataLen ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 20 19 22 0 ./CppStatUtilities.cc:140 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 54)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 22 20 23 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 23 22 24 1 ./CppStatUtilities.cc:140 (set (reg:DI 60 [ ivtmp.221 ])
        (reg/v/f:DI 63 [ data ])) -1 (nil)
    (nil))

(insn 24 23 25 1 ./CppStatUtilities.cc:140 (set (reg/v:SI 61 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 25 24 26 2 13 "" [1 uses])

(note 26 25 27 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 27 26 28 2 ("./CppStatUtilities.cc") 142)

(insn 28 27 29 2 ./CppStatUtilities.cc:142 (set (reg:DF 62 [ D.34589 ])
        (mem:DF (reg:DI 60 [ ivtmp.221 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 29 28 30 2 ./CppStatUtilities.cc:142 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 62 [ D.34589 ])
            (mem:DF (reg/v/f:DI 65 [ maxVal ]) [3 S8 A64]))) -1 (nil)
    (nil))

(jump_insn 30 29 56 2 ./CppStatUtilities.cc:142 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 34)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 56 30 31 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(jump_insn 31 56 32 3 ./CppStatUtilities.cc:142 (set (pc)
        (label_ref 40)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 32 31 34)

;; Start of basic block 4, registers live: (nil)
(code_label 34 32 35 4 16 "" [1 uses])

(note 35 34 36 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 36 35 37 4 ("./CppStatUtilities.cc") 144)

(insn 37 36 38 4 ./CppStatUtilities.cc:144 (set (mem:DF (reg/v/f:DI 65 [ maxVal ]) [3 S8 A64])
        (reg:DF 62 [ D.34589 ])) -1 (nil)
    (nil))

(note 38 37 39 4 ("./CppStatUtilities.cc") 145)

(insn 39 38 40 4 ./CppStatUtilities.cc:145 (set (mem:SI (reg/v/f:DI 66 [ maxInd ]) [5 S4 A32])
        (reg/v:SI 61 [ i ])) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 40 39 41 5 14 "" [1 uses])

(note 41 40 42 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 42 41 43 5 ("./CppStatUtilities.cc") 140)

(insn 43 42 44 5 ./CppStatUtilities.cc:140 (parallel [
            (set (reg/v:SI 61 [ i ])
                (plus:SI (reg/v:SI 61 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 44 43 45 5 ./CppStatUtilities.cc:140 (parallel [
            (set (reg:DI 60 [ ivtmp.221 ])
                (plus:DI (reg:DI 60 [ ivtmp.221 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 45 44 46 5 ./CppStatUtilities.cc:140 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 61 [ i ])
            (reg/v:SI 64 [ dataLen ]))) -1 (nil)
    (nil))

(jump_insn 46 45 49 5 ./CppStatUtilities.cc:140 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 25)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(note 49 46 50 NOTE_INSN_FUNCTION_END)

(note 50 49 54 ("./CppStatUtilities.cc") 148)

;; Start of basic block 6, registers live: (nil)
(code_label 54 50 59 6 17 "" [1 uses])

(note 59 54 0 6 [bb 6] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 6, registers live:
 (nil)


;; Function void whichMin(double*, int, double&, int&) (_Z8whichMinPdiRdRi)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->8 redirected to 10
Deleted label in block 2.
Edge 3->5 redirected to 6
Forwarding edge 3->4 to 7 failed.
Forwarding edge 3->4 to 7 failed.
Deleting block 5.
Forwarding edge 7->8 to 10 failed.
Deleted label in block 8.
Deleting fallthru block 8.
Forwarding edge 7->9 to 10 failed.
Deleting fallthru block 9.


try_optimize_cfg iteration 2

Forwarding edge 3->4 to 7 failed.


try_optimize_cfg iteration 1

Forwarding edge 2->3 to 5 failed.
(note 1 0 9 ("./CppStatUtilities.cc") 153)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 ./CppStatUtilities.cc:153 (set (reg/v/f:DI 61 [ data ])
        (reg:DI 5 di [ data ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:153 (set (reg/v:SI 62 [ dataLen ])
        (reg:SI 4 si [ dataLen ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:153 (set (reg/v/f:DI 63 [ minVal ])
        (reg:DI 1 dx [ minVal ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:153 (set (reg/v/f:DI 64 [ minInd ])
        (reg:DI 2 cx [ minInd ])) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("./CppStatUtilities.cc") 155)

(insn 12 11 13 0 ./CppStatUtilities.cc:155 (set (reg:DF 65)
        (mem:DF (reg/v/f:DI 61 [ data ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:155 (set (mem:DF (reg/v/f:DI 63 [ minVal ]) [3 S8 A64])
        (reg:DF 65)) -1 (nil)
    (nil))

(note 14 13 15 0 ("./CppStatUtilities.cc") 156)

(insn 15 14 16 0 ./CppStatUtilities.cc:156 (set (mem:SI (reg/v/f:DI 64 [ minInd ]) [5 S4 A32])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 16 15 17 0 ("./CppStatUtilities.cc") 159)

(insn 17 16 18 0 ./CppStatUtilities.cc:159 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 62 [ dataLen ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 18 17 20 0 ./CppStatUtilities.cc:159 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 53)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 20 18 21 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 21 20 22 1 ./CppStatUtilities.cc:159 (set (reg:DI 58 [ ivtmp.257 ])
        (reg/v/f:DI 61 [ data ])) -1 (nil)
    (nil))

(insn 22 21 23 1 ./CppStatUtilities.cc:159 (set (reg/v:SI 59 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 23 22 24 2 82 "" [1 uses])

(note 24 23 25 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 25 24 26 2 ("./CppStatUtilities.cc") 161)

(insn 26 25 27 2 ./CppStatUtilities.cc:161 (set (reg:DF 60 [ D.34607 ])
        (mem:DF (reg:DI 58 [ ivtmp.257 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 27 26 28 2 ./CppStatUtilities.cc:161 (set (reg:DF 66)
        (mem:DF (reg/v/f:DI 63 [ minVal ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 28 27 29 2 ./CppStatUtilities.cc:161 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 66)
            (reg:DF 60 [ D.34607 ]))) -1 (nil)
    (nil))

(jump_insn 29 28 55 2 ./CppStatUtilities.cc:161 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 33)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 55 29 30 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(jump_insn 30 55 31 3 ./CppStatUtilities.cc:161 (set (pc)
        (label_ref 39)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 31 30 33)

;; Start of basic block 4, registers live: (nil)
(code_label 33 31 34 4 85 "" [1 uses])

(note 34 33 35 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 35 34 36 4 ("./CppStatUtilities.cc") 163)

(insn 36 35 37 4 ./CppStatUtilities.cc:163 (set (mem:DF (reg/v/f:DI 63 [ minVal ]) [3 S8 A64])
        (reg:DF 60 [ D.34607 ])) -1 (nil)
    (nil))

(note 37 36 38 4 ("./CppStatUtilities.cc") 164)

(insn 38 37 39 4 ./CppStatUtilities.cc:164 (set (mem:SI (reg/v/f:DI 64 [ minInd ]) [5 S4 A32])
        (reg/v:SI 59 [ i ])) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 39 38 40 5 83 "" [1 uses])

(note 40 39 41 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 41 40 42 5 ("./CppStatUtilities.cc") 159)

(insn 42 41 43 5 ./CppStatUtilities.cc:159 (parallel [
            (set (reg/v:SI 59 [ i ])
                (plus:SI (reg/v:SI 59 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 43 42 44 5 ./CppStatUtilities.cc:159 (parallel [
            (set (reg:DI 58 [ ivtmp.257 ])
                (plus:DI (reg:DI 58 [ ivtmp.257 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 44 43 45 5 ./CppStatUtilities.cc:159 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 59 [ i ])
            (reg/v:SI 62 [ dataLen ]))) -1 (nil)
    (nil))

(jump_insn 45 44 48 5 ./CppStatUtilities.cc:159 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 23)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(note 48 45 49 NOTE_INSN_FUNCTION_END)

(note 49 48 53 ("./CppStatUtilities.cc") 167)

;; Start of basic block 6, registers live: (nil)
(code_label 53 49 58 6 86 "" [1 uses])

(note 58 53 0 6 [bb 6] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 6, registers live:
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

(note 1 0 7 ("./CppStatUtilities.cc") 601)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppStatUtilities.cc:601 (set (reg/v:SI 58 [ __initialize_p ])
        (reg:SI 5 di [ __initialize_p ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:601 (set (reg/v:SI 59 [ __priority ])
        (reg:SI 4 si [ __priority ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppStatUtilities.cc") 601)

(insn 10 9 11 0 ./CppStatUtilities.cc:601 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 58 [ __initialize_p ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./CppStatUtilities.cc:601 (set (pc)
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

(insn 14 13 15 1 ./CppStatUtilities.cc:601 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 59 [ __priority ])
            (const_int 65535 [0xffff]))) -1 (nil)
    (nil))

(jump_insn 15 14 17 1 ./CppStatUtilities.cc:601 (set (pc)
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
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2b5eb61c5c60 __ioinit>)) -1 (nil)
    (nil))

(call_insn 20 19 21 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitC1Ev") [flags 0x41] <function_decl 0x2b5eb598b000 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 1 dx)
        (symbol_ref:DI ("__dso_handle") [flags 0x40] <var_decl 0x2b5eb64fe840 __dso_handle>)) -1 (nil)
    (nil))

(insn 22 21 23 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("__tcf_0") [flags 0x3] <function_decl 0x2b5eb64c6400 __tcf_0>)) -1 (nil)
    (nil))

(call_insn/j 24 23 25 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_atexit") [flags 0x41] <function_decl 0x2b5eb64c6500 __cxa_atexit>) [0 S1 A8])
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

(note 29 28 33 ("./CppStatUtilities.cc") 601)

;; Start of basic block 3, registers live: (nil)
(code_label 33 29 36 3 152 "" [2 uses])

(note 36 33 0 3 [bb 3] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 3, registers live:
 (nil)


;; Function (static initializers for ./CppStatUtilities.cc) (_GLOBAL__I__ZN14summaryStats_tC2Ev)



try_optimize_cfg iteration 1

Deleting fallthru block 0.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 3 ("./CppStatUtilities.cc") 602)

(note 3 1 6 0 NOTE_INSN_FUNCTION_BEG)

;; Start of basic block 0, registers live: (nil)
(note 6 3 7 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(note 7 6 8 0 ("./CppStatUtilities.cc") 602)

(insn 8 7 9 0 ./CppStatUtilities.cc:602 (set (reg:SI 4 si)
        (const_int 65535 [0xffff])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppStatUtilities.cc:602 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 ./CppStatUtilities.cc:602 (call (mem:QI (symbol_ref:DI ("_Z41__static_initialization_and_destruction_0ii") [flags 0x3] <function_decl 0x2b5eb64c6300 __static_initialization_and_destruction_0>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("./CppStatUtilities.cc") 602)


;; Function void __tcf_0(void*) (__tcf_0)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg/v/f:DI 58 [ D.35343 ])
        (reg:DI 5 di [ D.35343 ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

(insn 9 8 10 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2b5eb61c5c60 __ioinit>)) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitD1Ev") [flags 0x41] <function_decl 0x2b5eb598b200 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)


;; Function void summaryStats_t::print(FILE*) const (_ZNK14summaryStats_t5printEP8_IO_FILE)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("./CppStatUtilities.cc") 63)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppStatUtilities.cc:63 (set (reg/f:DI 58 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:63 (set (reg/v/f:DI 59 [ f ])
        (reg:DI 4 si [ f ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppStatUtilities.cc") 65)

(insn 10 9 11 0 ./CppStatUtilities.cc:65 (set (reg:DF 60)
        (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppStatUtilities.cc:65 (set (reg:DF 61)
        (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppStatUtilities.cc:65 (set (reg:DF 62)
        (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:65 (set (reg:DF 63)
        (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppStatUtilities.cc:65 (set (reg:DF 64)
        (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:65 (set (reg:DF 65)
        (mem/s:DF (plus:DI (reg/f:DI 58 [ this ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppStatUtilities.cc:65 (set (reg:DF 27 xmm6)
        (reg:DF 60)) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:65 (set (reg:DF 26 xmm5)
        (reg:DF 61)) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:65 (set (reg:DF 25 xmm4)
        (reg:DF 62)) -1 (nil)
    (nil))

(insn 19 18 20 0 ./CppStatUtilities.cc:65 (set (reg:DF 24 xmm3)
        (reg:DF 63)) -1 (nil)
    (nil))

(insn 20 19 21 0 ./CppStatUtilities.cc:65 (set (reg:DF 23 xmm2)
        (reg:DF 64)) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppStatUtilities.cc:65 (set (reg:DF 22 xmm1)
        (reg:DF 65)) -1 (nil)
    (nil))

(insn 22 21 23 0 ./CppStatUtilities.cc:65 (set (reg:DF 21 xmm0)
        (mem/s:DF (reg/f:DI 58 [ this ]) [3 <variable>.min+0 S8 A64])) -1 (nil)
    (nil))

(insn 23 22 24 0 ./CppStatUtilities.cc:65 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC2") [flags 0x2] <string_cst 0x2b5eb6645740>)) -1 (nil)
    (nil))

(insn 24 23 25 0 ./CppStatUtilities.cc:65 (set (reg:DI 5 di)
        (reg/v/f:DI 59 [ f ])) -1 (nil)
    (nil))

(insn 25 24 26 0 ./CppStatUtilities.cc:65 (set (reg:QI 0 ax)
        (const_int 7 [0x7])) -1 (nil)
    (nil))

(call_insn/j 26 25 27 0 ./CppStatUtilities.cc:65 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
                    (expr_list:REG_DEP_TRUE (use (reg:DF 22 xmm1))
                        (expr_list:REG_DEP_TRUE (use (reg:DF 23 xmm2))
                            (expr_list:REG_DEP_TRUE (use (reg:DF 24 xmm3))
                                (expr_list:REG_DEP_TRUE (use (reg:DF 25 xmm4))
                                    (expr_list:REG_DEP_TRUE (use (reg:DF 26 xmm5))
                                        (expr_list:REG_DEP_TRUE (use (reg:DF 27 xmm6))
                                            (nil))))))))))))
;; End of basic block 0, registers live:
 (nil)

(barrier 27 26 28)

(note 28 27 29 NOTE_INSN_FUNCTION_END)

(note 29 28 0 ("./CppStatUtilities.cc") 66)


;; Function void signalPeakPar(const double*, int, PeakPar_t*&, int&) (_Z13signalPeakParPKdiRP9PeakPar_tRi)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Deleted label in block 6.
Deleted label in block 8.
Edge 8->10 redirected to 11
Forwarding edge 8->9 to 75 failed.
Forwarding edge 8->9 to 75 failed.
Deleting block 10.
Deleted label in block 13.
Edge 14->16 redirected to 17
Forwarding edge 14->15 to 18 failed.
Forwarding edge 14->15 to 18 failed.
Deleting block 16.
Deleted label in block 19.
Edge 20->22 redirected to 23
Forwarding edge 20->21 to 24 failed.
Forwarding edge 20->21 to 24 failed.
Deleting block 22.
Edge 24->26 redirected to 27
Forwarding edge 24->25 to 36 failed.
Forwarding edge 24->25 to 36 failed.
Deleting block 26.
Edge 27->29 redirected to 30
Forwarding edge 27->28 to 36 failed.
Forwarding edge 27->28 to 36 failed.
Deleting block 29.
Deleted label in block 31.
Deleted label in block 33.
Edge 36->38 redirected to 39
Forwarding edge 36->37 to 46 failed.
Forwarding edge 36->37 to 46 failed.
Deleting block 38.
Edge 39->41 redirected to 42
Forwarding edge 39->40 to 46 failed.
Forwarding edge 39->40 to 46 failed.
Deleting block 41.
Edge 46->48 redirected to 49
Forwarding edge 46->47 to 50 failed.
Forwarding edge 46->47 to 50 failed.
Deleting block 48.
Edge 50->52 redirected to 53
Forwarding edge 50->51 to 57 failed.
Forwarding edge 50->51 to 57 failed.
Deleting block 52.
Edge 53->55 redirected to 56
Forwarding edge 53->54 to 57 failed.
Forwarding edge 53->54 to 57 failed.
Deleting block 55.
Deleted label in block 58.
Edge 59->61 redirected to 62
Forwarding edge 60->61 to 62 failed.
Deleted label in block 61.
Deleting fallthru block 61.
Edge 63->65 redirected to 66
Forwarding edge 63->64 to 72 failed.
Forwarding edge 63->64 to 72 failed.
Deleting block 65.
Deleted label in block 69.
Forwarding edge 72->73 to 68 failed.


try_optimize_cfg iteration 2

Forwarding edge 8->9 to 75 failed.
Forwarding edge 14->15 to 18 failed.
Forwarding edge 20->21 to 24 failed.
Forwarding edge 24->25 to 36 failed.
Forwarding edge 27->28 to 36 failed.
Forwarding edge 36->37 to 46 failed.
Forwarding edge 39->40 to 46 failed.
Forwarding edge 46->47 to 50 failed.
Forwarding edge 50->51 to 57 failed.
Forwarding edge 53->54 to 57 failed.
Forwarding edge 63->64 to 72 failed.
Forwarding edge 72->73 to 68 failed.


try_optimize_cfg iteration 1

Forwarding edge 7->8 to 62 failed.
Forwarding edge 12->13 to 15 failed.
Forwarding edge 17->18 to 20 failed.
Forwarding edge 20->21 to 30 failed.
Forwarding edge 22->23 to 30 failed.
Forwarding edge 30->31 to 38 failed.
Forwarding edge 32->33 to 38 failed.
Forwarding edge 38->39 to 41 failed.
Forwarding edge 41->42 to 46 failed.
Forwarding edge 43->44 to 46 failed.
Forwarding edge 51->52 to 59 failed.
Forwarding edge 59->60 to 55 failed.
(note 1 0 26 ("./CppStatUtilities.cc") 498)

;; Start of basic block 0, registers live: (nil)
(note 26 1 20 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 20 26 21 0 ./CppStatUtilities.cc:498 (set (reg/v/f:DI 116 [ x ])
        (reg:DI 5 di [ x ])) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppStatUtilities.cc:498 (set (reg/v:SI 117 [ xSize ])
        (reg:SI 4 si [ xSize ])) -1 (nil)
    (nil))

(insn 22 21 23 0 ./CppStatUtilities.cc:498 (set (reg/v/f:DI 118 [ peakPar ])
        (reg:DI 1 dx [ peakPar ])) -1 (nil)
    (nil))

(insn 23 22 24 0 ./CppStatUtilities.cc:498 (set (reg/v/f:DI 119 [ nPeaks ])
        (reg:DI 2 cx [ nPeaks ])) -1 (nil)
    (nil))

(note 24 23 28 0 NOTE_INSN_FUNCTION_BEG)

(note 28 24 29 0 ("./CppStatUtilities.cc") 501)

(insn 29 28 30 0 ./CppStatUtilities.cc:501 (set (reg:DI 115 [ D.35195 ])
        (sign_extend:DI (reg/v:SI 117 [ xSize ]))) -1 (nil)
    (nil))

(insn 30 29 31 0 ./CppStatUtilities.cc:501 (set (reg:DI 4 si)
        (const_int 4 [0x4])) -1 (nil)
    (nil))

(insn 31 30 32 0 ./CppStatUtilities.cc:501 (set (reg:DI 5 di)
        (reg:DI 115 [ D.35195 ])) -1 (nil)
    (nil))

(call_insn 32 31 33 0 ./CppStatUtilities.cc:501 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b5eb3cbd200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 33 32 34 0 ./CppStatUtilities.cc:501 (set (reg/f:DI 120)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 120)
        (nil)))

(insn 34 33 35 0 ./CppStatUtilities.cc:501 (set (reg:DI 114 [ D.35196 ])
        (reg/f:DI 120)) -1 (nil)
    (nil))

(note 35 34 36 0 ("./CppStatUtilities.cc") 502)

(insn 36 35 37 0 ./CppStatUtilities.cc:502 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 114 [ D.35196 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 37 36 39 0 ./CppStatUtilities.cc:502 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 49)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 39 37 40 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 40 39 41 1 ./CppStatUtilities.cc:502 (set (reg:SI 2 cx)
        (const_int 502 [0x1f6])) -1 (nil)
    (nil))

(insn 41 40 42 1 ./CppStatUtilities.cc:502 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 42 41 43 1 ./CppStatUtilities.cc:502 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 43 42 44 1 ./CppStatUtilities.cc:502 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 44 43 45 1 ./CppStatUtilities.cc:502 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 45 44 46 1 ./CppStatUtilities.cc:502 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 46 45 47 1 ./CppStatUtilities.cc:502 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 47 46 48 1 ./CppStatUtilities.cc:502 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 1, registers live:
 (nil)

(barrier 48 47 49)

;; Start of basic block 2, registers live: (nil)
(code_label 49 48 50 2 161 "" [1 uses])

(note 50 49 51 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 51 50 52 2 ("./CppStatUtilities.cc") 504)

(insn 52 51 53 2 ./CppStatUtilities.cc:504 (set (reg:DI 4 si)
        (const_int 4 [0x4])) -1 (nil)
    (nil))

(insn 53 52 54 2 ./CppStatUtilities.cc:504 (set (reg:DI 5 di)
        (reg:DI 115 [ D.35195 ])) -1 (nil)
    (nil))

(call_insn 54 53 55 2 ./CppStatUtilities.cc:504 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b5eb3cbd200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 55 54 56 2 ./CppStatUtilities.cc:504 (set (reg/f:DI 121)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 121)
        (nil)))

(insn 56 55 57 2 ./CppStatUtilities.cc:504 (set (reg:DI 113 [ D.35200 ])
        (reg/f:DI 121)) -1 (nil)
    (nil))

(note 57 56 58 2 ("./CppStatUtilities.cc") 505)

(insn 58 57 59 2 ./CppStatUtilities.cc:505 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 113 [ D.35200 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 59 58 61 2 ./CppStatUtilities.cc:505 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 61 59 62 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 62 61 63 3 ./CppStatUtilities.cc:505 (set (reg:SI 2 cx)
        (const_int 505 [0x1f9])) -1 (nil)
    (nil))

(insn 63 62 64 3 ./CppStatUtilities.cc:505 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 64 63 65 3 ./CppStatUtilities.cc:505 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 65 64 66 3 ./CppStatUtilities.cc:505 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 66 65 67 3 ./CppStatUtilities.cc:505 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 67 66 68 3 ./CppStatUtilities.cc:505 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 68 67 69 3 ./CppStatUtilities.cc:505 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 69 68 70 3 ./CppStatUtilities.cc:505 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 3, registers live:
 (nil)

(barrier 70 69 71)

;; Start of basic block 4, registers live: (nil)
(code_label 71 70 72 4 163 "" [1 uses])

(note 72 71 73 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 73 72 74 4 ("./CppStatUtilities.cc") 507)

(insn 74 73 75 4 ./CppStatUtilities.cc:507 (set (reg:DI 4 si)
        (const_int 4 [0x4])) -1 (nil)
    (nil))

(insn 75 74 76 4 ./CppStatUtilities.cc:507 (set (reg:DI 5 di)
        (reg:DI 115 [ D.35195 ])) -1 (nil)
    (nil))

(call_insn 76 75 77 4 ./CppStatUtilities.cc:507 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b5eb3cbd200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 77 76 78 4 ./CppStatUtilities.cc:507 (set (reg/f:DI 122)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 122)
        (nil)))

(insn 78 77 79 4 ./CppStatUtilities.cc:507 (set (reg:DI 112 [ D.35203 ])
        (reg/f:DI 122)) -1 (nil)
    (nil))

(note 79 78 80 4 ("./CppStatUtilities.cc") 508)

(insn 80 79 81 4 ./CppStatUtilities.cc:508 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 112 [ D.35203 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 81 80 83 4 ./CppStatUtilities.cc:508 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 93)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 83 81 84 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 84 83 85 5 ./CppStatUtilities.cc:508 (set (reg:SI 2 cx)
        (const_int 508 [0x1fc])) -1 (nil)
    (nil))

(insn 85 84 86 5 ./CppStatUtilities.cc:508 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 86 85 87 5 ./CppStatUtilities.cc:508 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 87 86 88 5 ./CppStatUtilities.cc:508 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 88 87 89 5 ./CppStatUtilities.cc:508 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 89 88 90 5 ./CppStatUtilities.cc:508 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 90 89 91 5 ./CppStatUtilities.cc:508 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 91 90 92 5 ./CppStatUtilities.cc:508 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(barrier 92 91 93)

;; Start of basic block 6, registers live: (nil)
(code_label 93 92 94 6 165 "" [1 uses])

(note 94 93 95 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 95 94 96 6 ("./CppStatUtilities.cc") 501)

(insn 96 95 97 6 ./CppStatUtilities.cc:501 (set (reg/v/f:DI 107 [ minLInd ])
        (reg:DI 114 [ D.35196 ])) -1 (nil)
    (nil))

(note 97 96 98 6 ("./CppStatUtilities.cc") 504)

(insn 98 97 99 6 ./CppStatUtilities.cc:504 (set (reg/v/f:DI 106 [ maxInd ])
        (reg:DI 113 [ D.35200 ])) -1 (nil)
    (nil))

(note 99 98 100 6 ("./CppStatUtilities.cc") 507)

(insn 100 99 101 6 ./CppStatUtilities.cc:507 (set (reg/v/f:DI 105 [ minRInd ])
        (reg:DI 112 [ D.35203 ])) -1 (nil)
    (nil))

(note 101 100 102 6 ("./CppStatUtilities.cc") 512)

(insn 102 101 103 6 ./CppStatUtilities.cc:512 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 117 [ xSize ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 103 102 105 6 ./CppStatUtilities.cc:512 (set (pc)
        (if_then_else (le (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 449)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 2100 [0x834])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 105 103 106 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 106 105 107 7 ./CppStatUtilities.cc:512 (parallel [
            (set (reg:DI 111 [ D.35206 ])
                (plus:DI (reg/v/f:DI 116 [ x ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 107 106 108 7 ./CppStatUtilities.cc:512 (set (reg:DF 123)
        (mem:DF (reg:DI 111 [ D.35206 ]) [3 S8 A64])) -1 (nil)
    (nil))

(jump_insn 108 107 494 7 ./CppStatUtilities.cc:512 (parallel [
            (set (pc)
                (if_then_else (eq (reg:DF 123)
                        (mem:DF (reg/v/f:DI 116 [ x ]) [3 S8 A64]))
                    (label_ref:DI 112)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 494 108 109 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(jump_insn 109 494 110 8 ./CppStatUtilities.cc:512 (set (pc)
        (label_ref 464)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

(barrier 110 109 112)

;; Start of basic block 9, registers live: (nil)
(code_label 112 110 113 9 171 "" [1 uses])

(note 113 112 114 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 114 113 115 9 ./CppStatUtilities.cc:512 (set (reg:DI 86 [ ivtmp.495 ])
        (reg:DI 111 [ D.35206 ])) -1 (nil)
    (nil))

(insn 115 114 116 9 ./CppStatUtilities.cc:512 (set (reg/v:SI 104 [ i ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 116 115 117 10 172 "" [1 uses])

(note 117 116 118 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 118 117 119 10 ("./CppStatUtilities.cc") 514)

(insn 119 118 120 10 ./CppStatUtilities.cc:514 (parallel [
            (set (reg/v:SI 83 [ i.503 ])
                (plus:SI (reg/v:SI 104 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 120 119 121 10 ./CppStatUtilities.cc:514 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 83 [ i.503 ])
            (reg/v:SI 117 [ xSize ]))) -1 (nil)
    (nil))

(jump_insn 121 120 123 10 ./CppStatUtilities.cc:514 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 130)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 10, registers live:
 (nil)

;; Start of basic block 11, registers live: (nil)
(note 123 121 124 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 124 123 125 11 ./CppStatUtilities.cc:514 (set (reg/v:SI 82 [ i.504 ])
        (reg/v:SI 104 [ i ])) -1 (nil)
    (nil))

(insn 125 124 126 11 ./CppStatUtilities.cc:514 (set (reg:DI 124)
        (sign_extend:DI (reg/v:SI 82 [ i.504 ]))) -1 (nil)
    (nil))

(insn 126 125 127 11 ./CppStatUtilities.cc:514 (parallel [
            (set (reg:DI 125)
                (ashift:DI (reg:DI 124)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 127 126 128 11 ./CppStatUtilities.cc:514 (parallel [
            (set (reg:DI 93 [ prephitmp.466 ])
                (plus:DI (reg/v/f:DI 116 [ x ])
                    (reg:DI 125)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 128 127 129 11 ./CppStatUtilities.cc:514 (set (pc)
        (label_ref 146)) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

(barrier 129 128 130)

;; Start of basic block 12, registers live: (nil)
(code_label 130 129 131 12 173 "" [1 uses])

(note 131 130 132 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(note 132 131 133 12 ("./CppStatUtilities.cc") 498)

(insn 133 132 134 12 ./CppStatUtilities.cc:498 (set (reg:DI 93 [ prephitmp.466 ])
        (reg:DI 86 [ ivtmp.495 ])) -1 (nil)
    (nil))

(insn 134 133 135 12 ./CppStatUtilities.cc:498 (parallel [
            (set (reg:DI 86 [ ivtmp.495 ])
                (plus:DI (reg:DI 86 [ ivtmp.495 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 135 134 136 12 ("./CppStatUtilities.cc") 514)

(insn 136 135 137 12 ./CppStatUtilities.cc:514 (set (reg:DF 126)
        (mem:DF (plus:DI (reg:DI 93 [ prephitmp.466 ])
                (const_int 8 [0x8])) [3 S8 A64])) -1 (nil)
    (nil))

(jump_insn 137 136 496 12 ./CppStatUtilities.cc:514 (parallel [
            (set (pc)
                (if_then_else (eq (reg:DF 126)
                        (mem:DF (reg:DI 93 [ prephitmp.466 ]) [3 S8 A64]))
                    (label_ref:DI 141)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 8900 [0x22c4])
        (nil)))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(note 496 137 138 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(jump_insn 138 496 139 13 ./CppStatUtilities.cc:514 (set (pc)
        (label_ref 146)) -1 (nil)
    (nil))
;; End of basic block 13, registers live:
 (nil)

(barrier 139 138 141)

;; Start of basic block 14, registers live: (nil)
(code_label 141 139 142 14 177 "" [1 uses])

(note 142 141 143 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(insn 143 142 144 14 ./CppStatUtilities.cc:514 (set (reg/v:SI 104 [ i ])
        (reg/v:SI 83 [ i.503 ])) -1 (nil)
    (nil))

(jump_insn 144 143 145 14 ./CppStatUtilities.cc:514 (set (pc)
        (label_ref 116)) -1 (nil)
    (nil))
;; End of basic block 14, registers live:
 (nil)

(barrier 145 144 146)

;; Start of basic block 15, registers live: (nil)
(code_label 146 145 147 15 175 "" [2 uses])

(note 147 146 148 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 148 147 149 15 ("./CppStatUtilities.cc") 519)

(insn 149 148 150 15 ./CppStatUtilities.cc:519 (set (mem:SI (reg/v/f:DI 107 [ minLInd ]) [5 S4 A32])
        (reg/v:SI 104 [ i ])) -1 (nil)
    (nil))

(note 150 149 151 15 ("./CppStatUtilities.cc") 520)

(insn 151 150 152 15 ./CppStatUtilities.cc:520 (set (reg/v:DF 102 [ minVal ])
        (mem:DF (reg:DI 93 [ prephitmp.466 ]) [3 S8 A64])) -1 (nil)
    (nil))

(note 152 151 153 15 ("./CppStatUtilities.cc") 523)

(insn 153 152 154 15 ./CppStatUtilities.cc:523 (parallel [
            (set (reg/v:SI 85 [ i.501 ])
                (plus:SI (reg/v:SI 104 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 154 153 155 15 ("./CppStatUtilities.cc") 525)

(insn 155 154 156 15 ./CppStatUtilities.cc:525 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 117 [ xSize ])
            (reg/v:SI 85 [ i.501 ]))) -1 (nil)
    (nil))

(jump_insn 156 155 158 15 ./CppStatUtilities.cc:525 (set (pc)
        (if_then_else (le (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 477)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 0 [0x0])
        (nil)))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(note 158 156 159 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(insn 159 158 160 16 ./CppStatUtilities.cc:525 (set (reg/v:DF 101 [ maxVal ])
        (reg/v:DF 102 [ minVal ])) -1 (nil)
    (nil))

(jump_insn 160 159 161 16 ./CppStatUtilities.cc:525 (set (pc)
        (label_ref 472)) -1 (nil)
    (nil))
;; End of basic block 16, registers live:
 (nil)

(barrier 161 160 162)

;; Start of basic block 17, registers live: (nil)
(code_label 162 161 163 17 181 "" [2 uses])

(note 163 162 164 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(note 164 163 165 17 ("./CppStatUtilities.cc") 528)

(insn 165 164 166 17 ./CppStatUtilities.cc:528 (set (reg:DI 127)
        (sign_extend:DI (reg/v:SI 85 [ i.501 ]))) -1 (nil)
    (nil))

(insn 166 165 167 17 ./CppStatUtilities.cc:528 (parallel [
            (set (reg:DI 75 [ temp.534 ])
                (ashift:DI (reg:DI 127)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 167 166 168 17 ./CppStatUtilities.cc:528 (set (reg:DI 95 [ prephitmp.451 ])
        (reg:DI 75 [ temp.534 ])) -1 (nil)
    (nil))

(insn 168 167 169 17 ./CppStatUtilities.cc:528 (parallel [
            (set (reg:DI 94 [ prephitmp.453 ])
                (plus:DI (reg:DI 95 [ prephitmp.451 ])
                    (reg/v/f:DI 116 [ x ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 169 168 170 17 ./CppStatUtilities.cc:528 (set (reg/v:DF 76 [ minVal.533 ])
        (mem:DF (reg:DI 94 [ prephitmp.453 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 170 169 171 17 ./CppStatUtilities.cc:528 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 102 [ minVal ])
            (reg/v:DF 76 [ minVal.533 ]))) -1 (nil)
    (nil))

(jump_insn 171 170 498 17 ./CppStatUtilities.cc:528 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 175)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(note 498 171 172 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(jump_insn 172 498 173 18 ./CppStatUtilities.cc:528 (set (pc)
        (label_ref 180)) -1 (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

(barrier 173 172 175)

;; Start of basic block 19, registers live: (nil)
(code_label 175 173 176 19 184 "" [1 uses])

(note 176 175 177 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(insn 177 176 178 19 ./CppStatUtilities.cc:528 (set (reg/v:DF 102 [ minVal ])
        (reg/v:DF 76 [ minVal.533 ])) -1 (nil)
    (nil))

(jump_insn 178 177 179 19 ./CppStatUtilities.cc:528 (set (pc)
        (label_ref 280)) -1 (nil)
    (nil))
;; End of basic block 19, registers live:
 (nil)

(barrier 179 178 180)

;; Start of basic block 20, registers live: (nil)
(code_label 180 179 181 20 182 "" [1 uses])

(note 181 180 182 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(note 182 181 183 20 ("./CppStatUtilities.cc") 533)

(insn 183 182 184 20 ./CppStatUtilities.cc:533 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 76 [ minVal.533 ])
            (reg/v:DF 102 [ minVal ]))) -1 (nil)
    (nil))

(jump_insn 184 183 500 20 ./CppStatUtilities.cc:533 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 188)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 20, registers live:
 (nil)

;; Start of basic block 21, registers live: (nil)
(note 500 184 185 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(jump_insn 185 500 186 21 ./CppStatUtilities.cc:533 (set (pc)
        (label_ref 234)) -1 (nil)
    (nil))
;; End of basic block 21, registers live:
 (nil)

(barrier 186 185 188)

;; Start of basic block 22, registers live: (nil)
(code_label 188 186 189 22 188 "" [1 uses])

(note 189 188 190 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(insn 190 189 191 22 ./CppStatUtilities.cc:533 (parallel [
            (set (reg:DI 128)
                (plus:DI (reg/v/f:DI 116 [ x ])
                    (const_int -8 [0xfffffffffffffff8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 191 190 192 22 ./CppStatUtilities.cc:533 (parallel [
            (set (reg:DI 78 [ temp.527 ])
                (plus:DI (reg:DI 95 [ prephitmp.451 ])
                    (reg:DI 128)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 192 191 502 22 ./CppStatUtilities.cc:533 (parallel [
            (set (pc)
                (if_then_else (eq (reg/v:DF 102 [ minVal ])
                        (mem:DF (reg:DI 78 [ temp.527 ]) [3 S8 A64]))
                    (label_ref:DI 196)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 22, registers live:
 (nil)

;; Start of basic block 23, registers live: (nil)
(note 502 192 193 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(jump_insn 193 502 194 23 ./CppStatUtilities.cc:533 (set (pc)
        (label_ref 234)) -1 (nil)
    (nil))
;; End of basic block 23, registers live:
 (nil)

(barrier 194 193 196)

;; Start of basic block 24, registers live: (nil)
(code_label 196 194 197 24 190 "" [1 uses])

(note 197 196 198 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(note 198 197 199 24 ("./CppStatUtilities.cc") 535)

(insn 199 198 200 24 ./CppStatUtilities.cc:535 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 85 [ i.501 ])
            (const_int 2 [0x2]))) -1 (nil)
    (nil))

(jump_insn 200 199 202 24 ./CppStatUtilities.cc:535 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 206)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(note 202 200 203 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(insn 203 202 204 25 ./CppStatUtilities.cc:535 (parallel [
            (set (reg:SI 96 [ prephitmp.445 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 204 203 205 25 ./CppStatUtilities.cc:535 (set (pc)
        (label_ref 224)) -1 (nil)
    (nil))
;; End of basic block 25, registers live:
 (nil)

(barrier 205 204 206)

;; Start of basic block 26, registers live: (nil)
(code_label 206 205 207 26 191 "" [1 uses])

(note 207 206 208 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(insn 208 207 209 26 ./CppStatUtilities.cc:535 (set (reg:DF 129)
        (mem:DF (plus:DI (plus:DI (reg/v/f:DI 116 [ x ])
                    (reg:DI 95 [ prephitmp.451 ]))
                (const_int -16 [0xfffffffffffffff0])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 209 208 210 26 ./CppStatUtilities.cc:535 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 129)
            (reg/v:DF 102 [ minVal ]))) -1 (nil)
    (nil))

(jump_insn 210 209 212 26 ./CppStatUtilities.cc:535 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 216)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 26, registers live:
 (nil)

;; Start of basic block 27, registers live: (nil)
(note 212 210 213 27 [bb 27] NOTE_INSN_BASIC_BLOCK)

(insn 213 212 214 27 ./CppStatUtilities.cc:535 (parallel [
            (set (reg:SI 96 [ prephitmp.445 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 214 213 215 27 ./CppStatUtilities.cc:535 (set (pc)
        (label_ref 224)) -1 (nil)
    (nil))
;; End of basic block 27, registers live:
 (nil)

(barrier 215 214 216)

;; Start of basic block 28, registers live: (nil)
(code_label 216 215 217 28 194 "" [1 uses])

(note 217 216 218 28 [bb 28] NOTE_INSN_BASIC_BLOCK)

(note 218 217 219 28 ("./CppStatUtilities.cc") 537)

(insn 219 218 220 28 ./CppStatUtilities.cc:537 (parallel [
            (set (reg:SI 96 [ prephitmp.445 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 220 219 221 28 ./CppStatUtilities.cc:537 (set (reg:DI 130)
        (sign_extend:DI (reg/v:SI 103 [ count ]))) -1 (nil)
    (nil))

(insn 221 220 222 28 ./CppStatUtilities.cc:537 (set (mem:SI (plus:DI (mult:DI (reg:DI 130)
                    (const_int 4 [0x4]))
                (reg/v/f:DI 105 [ minRInd ])) [5 S4 A32])
        (reg:SI 96 [ prephitmp.445 ])) -1 (nil)
    (nil))

(note 222 221 223 28 ("./CppStatUtilities.cc") 538)

(insn 223 222 224 28 ./CppStatUtilities.cc:538 (parallel [
            (set (reg/v:SI 103 [ count ])
                (plus:SI (reg/v:SI 103 [ count ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 28, registers live:
 (nil)

;; Start of basic block 29, registers live: (nil)
(code_label 224 223 225 29 193 "" [2 uses])

(note 225 224 226 29 [bb 29] NOTE_INSN_BASIC_BLOCK)

(note 226 225 227 29 ("./CppStatUtilities.cc") 541)

(insn 227 226 228 29 ./CppStatUtilities.cc:541 (set (reg:DI 131)
        (sign_extend:DI (reg/v:SI 103 [ count ]))) -1 (nil)
    (nil))

(insn 228 227 229 29 ./CppStatUtilities.cc:541 (set (mem:SI (plus:DI (mult:DI (reg:DI 131)
                    (const_int 4 [0x4]))
                (reg/v/f:DI 107 [ minLInd ])) [5 S4 A32])
        (reg:SI 96 [ prephitmp.445 ])) -1 (nil)
    (nil))

(note 229 228 230 29 ("./CppStatUtilities.cc") 542)

(insn 230 229 231 29 ./CppStatUtilities.cc:542 (set (reg/v:DF 101 [ maxVal ])
        (mem:DF (reg:DI 78 [ temp.527 ]) [3 S8 A64])) -1 (nil)
    (nil))

(note 231 230 232 29 ("./CppStatUtilities.cc") 533)

(jump_insn 232 231 233 29 ./CppStatUtilities.cc:533 (set (pc)
        (label_ref 280)) -1 (nil)
    (nil))
;; End of basic block 29, registers live:
 (nil)

(barrier 233 232 234)

;; Start of basic block 30, registers live: (nil)
(code_label 234 233 235 30 186 "" [2 uses])

(note 235 234 236 30 [bb 30] NOTE_INSN_BASIC_BLOCK)

(note 236 235 237 30 ("./CppStatUtilities.cc") 545)

(jump_insn 237 236 504 30 ./CppStatUtilities.cc:545 (parallel [
            (set (pc)
                (if_then_else (eq (reg/v:DF 102 [ minVal ])
                        (reg/v:DF 76 [ minVal.533 ]))
                    (label_ref:DI 241)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 30, registers live:
 (nil)

;; Start of basic block 31, registers live: (nil)
(note 504 237 238 31 [bb 31] NOTE_INSN_BASIC_BLOCK)

(jump_insn 238 504 239 31 ./CppStatUtilities.cc:545 (set (pc)
        (label_ref 280)) -1 (nil)
    (nil))
;; End of basic block 31, registers live:
 (nil)

(barrier 239 238 241)

;; Start of basic block 32, registers live: (nil)
(code_label 241 239 242 32 197 "" [1 uses])

(note 242 241 243 32 [bb 32] NOTE_INSN_BASIC_BLOCK)

(insn 243 242 244 32 ./CppStatUtilities.cc:545 (parallel [
            (set (reg:DI 132)
                (plus:DI (reg/v/f:DI 116 [ x ])
                    (const_int -8 [0xfffffffffffffff8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 244 243 245 32 ./CppStatUtilities.cc:545 (parallel [
            (set (reg:DI 79 [ temp.526 ])
                (plus:DI (reg:DI 95 [ prephitmp.451 ])
                    (reg:DI 132)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 245 244 506 32 ./CppStatUtilities.cc:545 (parallel [
            (set (pc)
                (if_then_else (eq (reg/v:DF 102 [ minVal ])
                        (mem:DF (reg:DI 79 [ temp.526 ]) [3 S8 A64]))
                    (label_ref:DI 249)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 32, registers live:
 (nil)

;; Start of basic block 33, registers live: (nil)
(note 506 245 246 33 [bb 33] NOTE_INSN_BASIC_BLOCK)

(jump_insn 246 506 247 33 ./CppStatUtilities.cc:545 (set (pc)
        (label_ref 280)) -1 (nil)
    (nil))
;; End of basic block 33, registers live:
 (nil)

(barrier 247 246 249)

;; Start of basic block 34, registers live: (nil)
(code_label 249 247 250 34 199 "" [1 uses])

(note 250 249 251 34 [bb 34] NOTE_INSN_BASIC_BLOCK)

(note 251 250 252 34 ("./CppStatUtilities.cc") 547)

(insn 252 251 253 34 ./CppStatUtilities.cc:547 (set (reg:DI 133)
        (sign_extend:DI (reg/v:SI 103 [ count ]))) -1 (nil)
    (nil))

(insn 253 252 254 34 ./CppStatUtilities.cc:547 (parallel [
            (set (reg:SI 134)
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 254 253 255 34 ./CppStatUtilities.cc:547 (set (mem:SI (plus:DI (mult:DI (reg:DI 133)
                    (const_int 4 [0x4]))
                (reg/v/f:DI 105 [ minRInd ])) [5 S4 A32])
        (reg:SI 134)) -1 (nil)
    (expr_list:REG_EQUAL (plus:SI (reg/v:SI 85 [ i.501 ])
            (const_int -1 [0xffffffffffffffff]))
        (nil)))

(note 255 254 256 34 ("./CppStatUtilities.cc") 549)

(insn 256 255 257 34 ./CppStatUtilities.cc:549 (set (reg/v:DF 101 [ maxVal ])
        (mem:DF (reg:DI 79 [ temp.526 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 257 256 258 34 ./CppStatUtilities.cc:549 (parallel [
            (set (reg:DI 87 [ ivtmp.489 ])
                (plus:DI (reg:DI 75 [ temp.534 ])
                    (reg/v/f:DI 116 [ x ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 258 257 259 34 ./CppStatUtilities.cc:549 (set (pc)
        (label_ref 268)) -1 (nil)
    (nil))
;; End of basic block 34, registers live:
 (nil)

(barrier 259 258 260)

;; Start of basic block 35, registers live: (nil)
(code_label 260 259 261 35 201 "" [1 uses])

(note 261 260 262 35 [bb 35] NOTE_INSN_BASIC_BLOCK)

(note 262 261 263 35 ("./CppStatUtilities.cc") 553)

(insn 263 262 264 35 ./CppStatUtilities.cc:553 (parallel [
            (set (reg/v:SI 85 [ i.501 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 264 263 265 35 ./CppStatUtilities.cc:553 (parallel [
            (set (reg:DI 87 [ ivtmp.489 ])
                (plus:DI (reg:DI 87 [ ivtmp.489 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 265 264 266 35 ("./CppStatUtilities.cc") 551)

(insn 266 265 267 35 ./CppStatUtilities.cc:551 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 117 [ xSize ])
            (reg/v:SI 85 [ i.501 ]))) -1 (nil)
    (nil))

(jump_insn 267 266 268 35 ./CppStatUtilities.cc:551 (set (pc)
        (if_then_else (le (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 271)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1100 [0x44c])
        (nil)))
;; End of basic block 35, registers live:
 (nil)

;; Start of basic block 36, registers live: (nil)
(code_label 268 267 269 36 200 "" [1 uses])

(note 269 268 270 36 [bb 36] NOTE_INSN_BASIC_BLOCK)

(jump_insn 270 269 271 36 ./CppStatUtilities.cc:551 (parallel [
            (set (pc)
                (if_then_else (eq (reg/v:DF 102 [ minVal ])
                        (mem:DF (reg:DI 87 [ ivtmp.489 ]) [3 S8 A64]))
                    (label_ref 260)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9500 [0x251c])
        (nil)))
;; End of basic block 36, registers live:
 (nil)

;; Start of basic block 37, registers live: (nil)
(code_label 271 270 272 37 202 "" [1 uses])

(note 272 271 273 37 [bb 37] NOTE_INSN_BASIC_BLOCK)

(note 273 272 274 37 ("./CppStatUtilities.cc") 548)

(insn 274 273 275 37 ./CppStatUtilities.cc:548 (parallel [
            (set (reg/v:SI 103 [ count ])
                (plus:SI (reg/v:SI 103 [ count ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 275 274 276 37 ("./CppStatUtilities.cc") 555)

(insn 276 275 277 37 ./CppStatUtilities.cc:555 (parallel [
            (set (reg/v:SI 85 [ i.501 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 277 276 278 37 ./CppStatUtilities.cc:555 (set (reg:DI 135)
        (sign_extend:DI (reg/v:SI 85 [ i.501 ]))) -1 (nil)
    (nil))

(insn 278 277 279 37 ./CppStatUtilities.cc:555 (parallel [
            (set (reg:DI 95 [ prephitmp.451 ])
                (ashift:DI (reg:DI 135)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 279 278 280 37 ./CppStatUtilities.cc:555 (parallel [
            (set (reg:DI 94 [ prephitmp.453 ])
                (plus:DI (reg/v/f:DI 116 [ x ])
                    (reg:DI 95 [ prephitmp.451 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 37, registers live:
 (nil)

;; Start of basic block 38, registers live: (nil)
(code_label 280 279 281 38 185 "" [4 uses])

(note 281 280 282 38 [bb 38] NOTE_INSN_BASIC_BLOCK)

(note 282 281 283 38 ("./CppStatUtilities.cc") 559)

(insn 283 282 284 38 ./CppStatUtilities.cc:559 (set (reg/v:DF 77 [ maxVal.532 ])
        (mem:DF (reg:DI 94 [ prephitmp.453 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 284 283 285 38 ./CppStatUtilities.cc:559 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 77 [ maxVal.532 ])
            (reg/v:DF 101 [ maxVal ]))) -1 (nil)
    (nil))

(jump_insn 285 284 508 38 ./CppStatUtilities.cc:559 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 289)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 38, registers live:
 (nil)

;; Start of basic block 39, registers live: (nil)
(note 508 285 286 39 [bb 39] NOTE_INSN_BASIC_BLOCK)

(jump_insn 286 508 287 39 ./CppStatUtilities.cc:559 (set (pc)
        (label_ref 294)) -1 (nil)
    (nil))
;; End of basic block 39, registers live:
 (nil)

(barrier 287 286 289)

;; Start of basic block 40, registers live: (nil)
(code_label 289 287 290 40 205 "" [1 uses])

(note 290 289 291 40 [bb 40] NOTE_INSN_BASIC_BLOCK)

(insn 291 290 292 40 ./CppStatUtilities.cc:559 (set (reg/v:DF 101 [ maxVal ])
        (reg/v:DF 77 [ maxVal.532 ])) -1 (nil)
    (nil))

(jump_insn 292 291 293 40 ./CppStatUtilities.cc:559 (set (pc)
        (label_ref 316)) -1 (nil)
    (nil))
;; End of basic block 40, registers live:
 (nil)

(barrier 293 292 294)

;; Start of basic block 41, registers live: (nil)
(code_label 294 293 295 41 203 "" [1 uses])

(note 295 294 296 41 [bb 41] NOTE_INSN_BASIC_BLOCK)

(note 296 295 297 41 ("./CppStatUtilities.cc") 563)

(insn 297 296 298 41 ./CppStatUtilities.cc:563 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 101 [ maxVal ])
            (reg/v:DF 77 [ maxVal.532 ]))) -1 (nil)
    (nil))

(jump_insn 298 297 510 41 ./CppStatUtilities.cc:563 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 302)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 41, registers live:
 (nil)

;; Start of basic block 42, registers live: (nil)
(note 510 298 299 42 [bb 42] NOTE_INSN_BASIC_BLOCK)

(jump_insn 299 510 300 42 ./CppStatUtilities.cc:563 (set (pc)
        (label_ref 316)) -1 (nil)
    (nil))
;; End of basic block 42, registers live:
 (nil)

(barrier 300 299 302)

;; Start of basic block 43, registers live: (nil)
(code_label 302 300 303 43 208 "" [1 uses])

(note 303 302 304 43 [bb 43] NOTE_INSN_BASIC_BLOCK)

(jump_insn 304 303 512 43 ./CppStatUtilities.cc:563 (parallel [
            (set (pc)
                (if_then_else (eq (reg/v:DF 101 [ maxVal ])
                        (mem:DF (plus:DI (plus:DI (reg/v/f:DI 116 [ x ])
                                    (reg:DI 95 [ prephitmp.451 ]))
                                (const_int -8 [0xfffffffffffffff8])) [3 S8 A64]))
                    (label_ref:DI 308)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 43, registers live:
 (nil)

;; Start of basic block 44, registers live: (nil)
(note 512 304 305 44 [bb 44] NOTE_INSN_BASIC_BLOCK)

(jump_insn 305 512 306 44 ./CppStatUtilities.cc:563 (set (pc)
        (label_ref 316)) -1 (nil)
    (nil))
;; End of basic block 44, registers live:
 (nil)

(barrier 306 305 308)

;; Start of basic block 45, registers live: (nil)
(code_label 308 306 309 45 210 "" [1 uses])

(note 309 308 310 45 [bb 45] NOTE_INSN_BASIC_BLOCK)

(note 310 309 311 45 ("./CppStatUtilities.cc") 565)

(insn 311 310 312 45 ./CppStatUtilities.cc:565 (set (reg:DI 136)
        (sign_extend:DI (reg/v:SI 103 [ count ]))) -1 (nil)
    (nil))

(insn 312 311 313 45 ./CppStatUtilities.cc:565 (parallel [
            (set (reg:SI 137)
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 313 312 314 45 ./CppStatUtilities.cc:565 (set (mem:SI (plus:DI (mult:DI (reg:DI 136)
                    (const_int 4 [0x4]))
                (reg/v/f:DI 106 [ maxInd ])) [5 S4 A32])
        (reg:SI 137)) -1 (nil)
    (expr_list:REG_EQUAL (plus:SI (reg/v:SI 85 [ i.501 ])
            (const_int -1 [0xffffffffffffffff]))
        (nil)))

(note 314 313 315 45 ("./CppStatUtilities.cc") 566)

(insn 315 314 316 45 ./CppStatUtilities.cc:566 (set (reg/v:DF 102 [ minVal ])
        (mem:DF (reg:DI 94 [ prephitmp.453 ]) [3 S8 A64])) -1 (nil)
    (nil))
;; End of basic block 45, registers live:
 (nil)

;; Start of basic block 46, registers live: (nil)
(code_label 316 315 317 46 206 "" [3 uses])

(note 317 316 318 46 [bb 46] NOTE_INSN_BASIC_BLOCK)

(note 318 317 319 46 ("./CppStatUtilities.cc") 568)

(insn 319 318 320 46 ./CppStatUtilities.cc:568 (parallel [
            (set (reg/v:SI 85 [ i.501 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 320 319 321 46 ("./CppStatUtilities.cc") 525)

(insn 321 320 322 46 ./CppStatUtilities.cc:525 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 117 [ xSize ])
            (reg/v:SI 85 [ i.501 ]))) -1 (nil)
    (nil))

(jump_insn 322 321 324 46 ./CppStatUtilities.cc:525 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 162)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 46, registers live:
 (nil)

;; Start of basic block 47, registers live: (nil)
(note 324 322 325 47 [bb 47] NOTE_INSN_BASIC_BLOCK)

(insn 325 324 326 47 ./CppStatUtilities.cc:525 (set (reg/v:SI 80 [ count.510 ])
        (reg/v:SI 103 [ count ])) -1 (nil)
    (nil))

(insn 326 325 327 47 ./CppStatUtilities.cc:525 (parallel [
            (set (reg/v:SI 84 [ i.502 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 327 326 328 47 ./CppStatUtilities.cc:525 (set (reg:DI 138)
        (sign_extend:DI (reg/v:SI 80 [ count.510 ]))) -1 (nil)
    (nil))

(insn 328 327 329 47 ./CppStatUtilities.cc:525 (parallel [
            (set (reg:DI 98 [ prephitmp.418 ])
                (ashift:DI (reg:DI 138)
                    (const_int 2 [0x2])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 329 328 330 47 ./CppStatUtilities.cc:525 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 80 [ count.510 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(insn 330 329 331 47 ./CppStatUtilities.cc:525 (set (reg:QI 139)
        (eq:QI (reg:CCZ 17 flags)
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EQUAL (eq:QI (reg/v:SI 80 [ count.510 ])
            (const_int 0 [0x0]))
        (nil)))

(insn 331 330 332 47 ./CppStatUtilities.cc:525 (set (reg:QI 97 [ prephitmp.420 ])
        (reg:QI 139)) -1 (nil)
    (nil))

(insn 332 331 333 47 ./CppStatUtilities.cc:525 (parallel [
            (set (reg:SI 100 [ prephitmp.406 ])
                (plus:SI (reg/v:SI 80 [ count.510 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 333 332 334 47 ./CppStatUtilities.cc:525 (set (reg:DI 140)
        (sign_extend:DI (reg:SI 100 [ prephitmp.406 ]))) -1 (nil)
    (nil))

(insn 334 333 335 47 ./CppStatUtilities.cc:525 (set (reg:DI 141)
        (reg:DI 140)) -1 (nil)
    (nil))

(insn 335 334 336 47 ./CppStatUtilities.cc:525 (parallel [
            (set (reg:DI 142)
                (ashift:DI (reg:DI 141)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 336 335 337 47 ./CppStatUtilities.cc:525 (parallel [
            (set (reg:DI 143)
                (plus:DI (reg:DI 142)
                    (reg:DI 140)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (mult:DI (reg:DI 140)
            (const_int 3 [0x3]))
        (nil)))

(insn 337 336 338 47 ./CppStatUtilities.cc:525 (parallel [
            (set (reg:DI 144)
                (ashift:DI (reg:DI 143)
                    (const_int 4 [0x4])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (mult:DI (reg:DI 140)
            (const_int 48 [0x30]))
        (nil)))

(insn 338 337 339 47 ./CppStatUtilities.cc:525 (set (reg:DI 99 [ prephitmp.410 ])
        (reg:DI 144)) -1 (nil)
    (nil))
;; End of basic block 47, registers live:
 (nil)

;; Start of basic block 48, registers live: (nil)
(code_label 339 338 340 48 212 "" [1 uses])

(note 340 339 341 48 [bb 48] NOTE_INSN_BASIC_BLOCK)

(note 341 340 342 48 ("./CppStatUtilities.cc") 573)

(insn 342 341 343 48 ./CppStatUtilities.cc:573 (set (reg:CCZ 17 flags)
        (compare:CCZ (mem:SI (plus:DI (reg/v/f:DI 107 [ minLInd ])
                    (reg:DI 98 [ prephitmp.418 ])) [5 S4 A32])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 343 342 514 48 ./CppStatUtilities.cc:573 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 347)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 48, registers live:
 (nil)

;; Start of basic block 49, registers live: (nil)
(note 514 343 344 49 [bb 49] NOTE_INSN_BASIC_BLOCK)

(insn 344 514 345 49 ./CppStatUtilities.cc:573 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:QI 97 [ prephitmp.420 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 345 344 347 49 ./CppStatUtilities.cc:573 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 378)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 49, registers live:
 (nil)

;; Start of basic block 50, registers live: (nil)
(code_label 347 345 348 50 215 "" [1 uses])

(note 348 347 349 50 [bb 50] NOTE_INSN_BASIC_BLOCK)

(insn 349 348 350 50 ./CppStatUtilities.cc:573 (parallel [
            (set (reg:SI 100 [ prephitmp.406 ])
                (plus:SI (reg/v:SI 103 [ count ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 350 349 351 50 ./CppStatUtilities.cc:573 (set (reg:DI 145)
        (sign_extend:DI (reg:SI 100 [ prephitmp.406 ]))) -1 (nil)
    (nil))

(insn 351 350 352 50 ./CppStatUtilities.cc:573 (set (reg:DI 146)
        (reg:DI 145)) -1 (nil)
    (nil))

(insn 352 351 353 50 ./CppStatUtilities.cc:573 (parallel [
            (set (reg:DI 147)
                (ashift:DI (reg:DI 146)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 353 352 354 50 ./CppStatUtilities.cc:573 (parallel [
            (set (reg:DI 148)
                (plus:DI (reg:DI 147)
                    (reg:DI 145)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (mult:DI (reg:DI 145)
            (const_int 3 [0x3]))
        (nil)))

(insn 354 353 355 50 ./CppStatUtilities.cc:573 (parallel [
            (set (reg:DI 149)
                (ashift:DI (reg:DI 148)
                    (const_int 4 [0x4])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (mult:DI (reg:DI 145)
            (const_int 48 [0x30]))
        (nil)))

(insn 355 354 356 50 ./CppStatUtilities.cc:573 (set (reg:DI 99 [ prephitmp.410 ])
        (reg:DI 149)) -1 (nil)
    (nil))
;; End of basic block 50, registers live:
 (nil)

;; Start of basic block 51, registers live: (nil)
(code_label 356 355 357 51 216 "" [1 uses])

(note 357 356 358 51 [bb 51] NOTE_INSN_BASIC_BLOCK)

(insn 358 357 359 51 ./CppStatUtilities.cc:573 (set (reg:DI 150)
        (sign_extend:DI (reg/v:SI 84 [ i.502 ]))) -1 (nil)
    (nil))

(insn 359 358 360 51 ./CppStatUtilities.cc:573 (set (reg:DF 110 [ D.35219 ])
        (mem:DF (plus:DI (mult:DI (reg:DI 150)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 116 [ x ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 360 359 361 51 ./CppStatUtilities.cc:573 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 110 [ D.35219 ])
            (reg/v:DF 101 [ maxVal ]))) -1 (nil)
    (nil))

(jump_insn 361 360 516 51 ./CppStatUtilities.cc:573 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 365)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 51, registers live:
 (nil)

;; Start of basic block 52, registers live: (nil)
(note 516 361 362 52 [bb 52] NOTE_INSN_BASIC_BLOCK)

(jump_insn 362 516 363 52 ./CppStatUtilities.cc:573 (set (pc)
        (label_ref 441)) -1 (nil)
    (nil))
;; End of basic block 52, registers live:
 (nil)

(barrier 363 362 365)

;; Start of basic block 53, registers live: (nil)
(code_label 365 363 366 53 219 "" [1 uses])

(note 366 365 367 53 [bb 53] NOTE_INSN_BASIC_BLOCK)

(note 367 366 368 53 ("./CppStatUtilities.cc") 575)

(insn 368 367 369 53 ./CppStatUtilities.cc:575 (set (mem:SI (plus:DI (reg/v/f:DI 106 [ maxInd ])
                (reg:DI 98 [ prephitmp.418 ])) [5 S4 A32])
        (reg/v:SI 84 [ i.502 ])) -1 (nil)
    (nil))

(note 369 368 370 53 ("./CppStatUtilities.cc") 576)

(insn 370 369 371 53 ./CppStatUtilities.cc:576 (set (mem:SI (plus:DI (reg/v/f:DI 105 [ minRInd ])
                (reg:DI 98 [ prephitmp.418 ])) [5 S4 A32])
        (reg/v:SI 84 [ i.502 ])) -1 (nil)
    (nil))

(note 371 370 372 53 ("./CppStatUtilities.cc") 573)

(jump_insn 372 371 373 53 ./CppStatUtilities.cc:573 (set (pc)
        (label_ref 378)) -1 (nil)
    (nil))
;; End of basic block 53, registers live:
 (nil)

(barrier 373 372 374)

;; Start of basic block 54, registers live: (nil)
(code_label 374 373 375 54 220 "" [1 uses])

(note 375 374 376 54 [bb 54] NOTE_INSN_BASIC_BLOCK)

(note 376 375 377 54 ("./CppStatUtilities.cc") 580)

(insn 377 376 378 54 ./CppStatUtilities.cc:580 (set (mem:SI (plus:DI (reg/v/f:DI 105 [ minRInd ])
                (reg:DI 98 [ prephitmp.418 ])) [5 S4 A32])
        (reg/v:SI 84 [ i.502 ])) -1 (nil)
    (nil))
;; End of basic block 54, registers live:
 (nil)

;; Start of basic block 55, registers live: (nil)
(code_label 378 377 379 55 213 "" [3 uses])

(note 379 378 380 55 [bb 55] NOTE_INSN_BASIC_BLOCK)

(note 380 379 381 55 ("./CppStatUtilities.cc") 583)

(insn 381 380 382 55 ./CppStatUtilities.cc:583 (set (mem:SI (reg/v/f:DI 119 [ nPeaks ]) [5 S4 A32])
        (reg:SI 100 [ prephitmp.406 ])) -1 (nil)
    (nil))

(note 382 381 383 55 ("./CppStatUtilities.cc") 584)

(insn 383 382 384 55 ./CppStatUtilities.cc:584 (set (reg:DI 5 di)
        (reg:DI 99 [ prephitmp.410 ])) -1 (nil)
    (nil))

(call_insn 384 383 385 55 ./CppStatUtilities.cc:584 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 385 384 386 55 ./CppStatUtilities.cc:584 (set (reg/f:DI 151)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 151)
        (nil)))

(insn 386 385 387 55 ./CppStatUtilities.cc:584 (set (reg:DI 91 [ ivtmp.477 ])
        (reg/f:DI 151)) -1 (nil)
    (nil))

(insn 387 386 388 55 ./CppStatUtilities.cc:584 (set (mem/f:DI (reg/v/f:DI 118 [ peakPar ]) [17 S8 A64])
        (reg:DI 91 [ ivtmp.477 ])) -1 (nil)
    (nil))

(note 388 387 389 55 ("./CppStatUtilities.cc") 586)

(insn 389 388 390 55 ./CppStatUtilities.cc:586 (set (reg:CCNO 17 flags)
        (compare:CCNO (mem:SI (reg/v/f:DI 119 [ nPeaks ]) [5 S4 A32])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 390 389 392 55 ./CppStatUtilities.cc:586 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 429)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 55, registers live:
 (nil)

;; Start of basic block 56, registers live: (nil)
(note 392 390 393 56 [bb 56] NOTE_INSN_BASIC_BLOCK)

(insn 393 392 394 56 ./CppStatUtilities.cc:586 (set (reg:DI 92 [ ivtmp.471 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 56, registers live:
 (nil)

;; Start of basic block 57, registers live: (nil)
(code_label 394 393 395 57 223 "" [1 uses])

(note 395 394 396 57 [bb 57] NOTE_INSN_BASIC_BLOCK)

(note 396 395 397 57 ("./CppStatUtilities.cc") 498)

(insn 397 396 398 57 ./CppStatUtilities.cc:498 (set (reg:DI 90 [ D.36135 ])
        (reg:DI 92 [ ivtmp.471 ])) -1 (nil)
    (nil))

(insn 398 397 399 57 ./CppStatUtilities.cc:498 (set (reg:DI 89 [ D.36136 ])
        (reg:DI 91 [ ivtmp.477 ])) -1 (nil)
    (nil))

(note 399 398 400 57 ("./CppStatUtilities.cc") 588)

(insn 400 399 401 57 ./CppStatUtilities.cc:588 (set (reg:SI 153)
        (mem:SI (plus:DI (mult:DI (reg:DI 90 [ D.36135 ])
                    (const_int 4 [0x4]))
                (reg/v/f:DI 107 [ minLInd ])) [5 S4 A32])) -1 (nil)
    (nil))

(insn 401 400 402 57 ./CppStatUtilities.cc:588 (set (reg:DI 152)
        (sign_extend:DI (reg:SI 153))) -1 (nil)
    (nil))

(insn 402 401 403 57 ./CppStatUtilities.cc:588 (set (reg:DF 154)
        (mem:DF (plus:DI (mult:DI (reg:DI 152)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 116 [ x ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 403 402 404 57 ./CppStatUtilities.cc:588 (set (mem/s:DF (plus:DI (reg:DI 89 [ D.36136 ])
                (const_int 24 [0x18])) [3 <variable>.startVal+0 S8 A8])
        (reg:DF 154)) -1 (nil)
    (nil))

(note 404 403 405 57 ("./CppStatUtilities.cc") 498)

(insn 405 404 406 57 ./CppStatUtilities.cc:498 (set (reg:DI 88 [ D.36139 ])
        (reg:DI 91 [ ivtmp.477 ])) -1 (nil)
    (nil))

(note 406 405 407 57 ("./CppStatUtilities.cc") 589)

(insn 407 406 408 57 ./CppStatUtilities.cc:589 (set (reg:SI 155)
        (mem:SI (plus:DI (mult:DI (reg:DI 90 [ D.36135 ])
                    (const_int 4 [0x4]))
                (reg/v/f:DI 107 [ minLInd ])) [5 S4 A32])) -1 (nil)
    (nil))

(insn 408 407 409 57 ./CppStatUtilities.cc:589 (set (mem/s:SI (plus:DI (reg:DI 88 [ D.36139 ])
                (const_int 16 [0x10])) [5 <variable>.startInd+0 S4 A8])
        (reg:SI 155)) -1 (nil)
    (nil))

(note 409 408 410 57 ("./CppStatUtilities.cc") 591)

(insn 410 409 411 57 ./CppStatUtilities.cc:591 (set (reg:SI 109 [ D.35278 ])
        (mem:SI (plus:DI (mult:DI (reg:DI 90 [ D.36135 ])
                    (const_int 4 [0x4]))
                (reg/v/f:DI 105 [ minRInd ])) [5 S4 A32])) -1 (nil)
    (nil))

(insn 411 410 412 57 ./CppStatUtilities.cc:591 (set (reg:DI 156)
        (sign_extend:DI (reg:SI 109 [ D.35278 ]))) -1 (nil)
    (nil))

(insn 412 411 413 57 ./CppStatUtilities.cc:591 (set (reg:DF 157)
        (mem:DF (plus:DI (mult:DI (reg:DI 156)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 116 [ x ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 413 412 414 57 ./CppStatUtilities.cc:591 (set (mem/s:DF (plus:DI (reg:DI 89 [ D.36136 ])
                (const_int 40 [0x28])) [3 <variable>.endVal+0 S8 A8])
        (reg:DF 157)) -1 (nil)
    (nil))

(note 414 413 415 57 ("./CppStatUtilities.cc") 592)

(insn 415 414 416 57 ./CppStatUtilities.cc:592 (set (mem/s:SI (plus:DI (reg:DI 88 [ D.36139 ])
                (const_int 32 [0x20])) [5 <variable>.endInd+0 S4 A8])
        (reg:SI 109 [ D.35278 ])) -1 (nil)
    (nil))

(note 416 415 417 57 ("./CppStatUtilities.cc") 594)

(insn 417 416 418 57 ./CppStatUtilities.cc:594 (set (reg:SI 108 [ D.35285 ])
        (mem:SI (plus:DI (mult:DI (reg:DI 90 [ D.36135 ])
                    (const_int 4 [0x4]))
                (reg/v/f:DI 106 [ maxInd ])) [5 S4 A32])) -1 (nil)
    (nil))

(insn 418 417 419 57 ./CppStatUtilities.cc:594 (set (reg:DI 158)
        (sign_extend:DI (reg:SI 108 [ D.35285 ]))) -1 (nil)
    (nil))

(insn 419 418 420 57 ./CppStatUtilities.cc:594 (set (reg:DF 159)
        (mem:DF (plus:DI (mult:DI (reg:DI 158)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 116 [ x ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 420 419 421 57 ./CppStatUtilities.cc:594 (set (mem/s:DF (plus:DI (reg:DI 89 [ D.36136 ])
                (const_int 8 [0x8])) [3 <variable>.maxVal+0 S8 A8])
        (reg:DF 159)) -1 (nil)
    (nil))

(note 421 420 422 57 ("./CppStatUtilities.cc") 595)

(insn 422 421 423 57 ./CppStatUtilities.cc:595 (set (mem/s:SI (reg:DI 88 [ D.36139 ]) [5 <variable>.maxInd+0 S4 A8])
        (reg:SI 108 [ D.35285 ])) -1 (nil)
    (nil))

(insn 423 422 424 57 ./CppStatUtilities.cc:595 (parallel [
            (set (reg/v:SI 81 [ i.505 ])
                (plus:SI (subreg:SI (reg:DI 92 [ ivtmp.471 ]) 0)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 424 423 425 57 ./CppStatUtilities.cc:595 (parallel [
            (set (reg:DI 92 [ ivtmp.471 ])
                (plus:DI (reg:DI 92 [ ivtmp.471 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 425 424 426 57 ./CppStatUtilities.cc:595 (parallel [
            (set (reg:DI 91 [ ivtmp.477 ])
                (plus:DI (reg:DI 91 [ ivtmp.477 ])
                    (const_int 48 [0x30])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 426 425 427 57 ("./CppStatUtilities.cc") 586)

(insn 427 426 428 57 ./CppStatUtilities.cc:586 (set (reg:CCGC 17 flags)
        (compare:CCGC (mem:SI (reg/v/f:DI 119 [ nPeaks ]) [5 S4 A32])
            (reg/v:SI 81 [ i.505 ]))) -1 (nil)
    (nil))

(jump_insn 428 427 429 57 ./CppStatUtilities.cc:586 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 394)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 57, registers live:
 (nil)

;; Start of basic block 58, registers live: (nil)
(code_label 429 428 430 58 221 "" [1 uses])

(note 430 429 431 58 [bb 58] NOTE_INSN_BASIC_BLOCK)

(note 431 430 432 58 ("./CppStatUtilities.cc") 598)

(insn 432 431 433 58 ./CppStatUtilities.cc:598 (set (reg:DI 5 di)
        (reg/v/f:DI 107 [ minLInd ])) -1 (nil)
    (nil))

(call_insn 433 432 434 58 ./CppStatUtilities.cc:598 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 434 433 435 58 ("./CppStatUtilities.cc") 599)

(insn 435 434 436 58 ./CppStatUtilities.cc:599 (set (reg:DI 5 di)
        (reg/v/f:DI 105 [ minRInd ])) -1 (nil)
    (nil))

(call_insn 436 435 437 58 ./CppStatUtilities.cc:599 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 437 436 438 58 ("./CppStatUtilities.cc") 600)

(insn 438 437 439 58 ./CppStatUtilities.cc:600 (set (reg:DI 5 di)
        (reg/v/f:DI 106 [ maxInd ])) -1 (nil)
    (nil))

(call_insn/j 439 438 440 58 ./CppStatUtilities.cc:600 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 58, registers live:
 (nil)

(barrier 440 439 441)

;; Start of basic block 59, registers live: (nil)
(code_label 441 440 442 59 217 "" [1 uses])

(note 442 441 443 59 [bb 59] NOTE_INSN_BASIC_BLOCK)

(note 443 442 444 59 ("./CppStatUtilities.cc") 578)

(insn 444 443 445 59 ./CppStatUtilities.cc:578 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 102 [ minVal ])
            (reg:DF 110 [ D.35219 ]))) -1 (nil)
    (nil))

(jump_insn 445 444 448 59 ./CppStatUtilities.cc:578 (set (pc)
        (if_then_else (ge (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 374)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 59, registers live:
 (nil)

;; Start of basic block 60, registers live: (nil)
(note 448 445 446 60 [bb 60] NOTE_INSN_BASIC_BLOCK)

(jump_insn 446 448 447 60 ./CppStatUtilities.cc:578 (set (pc)
        (label_ref 378)) -1 (nil)
    (nil))
;; End of basic block 60, registers live:
 (nil)

(barrier 447 446 449)

;; Start of basic block 61, registers live: (nil)
(code_label 449 447 450 61 167 "" [1 uses])

(note 450 449 451 61 [bb 61] NOTE_INSN_BASIC_BLOCK)

(note 451 450 452 61 ("./CppStatUtilities.cc") 519)

(insn 452 451 453 61 ./CppStatUtilities.cc:519 (set (mem:SI (reg/v/f:DI 107 [ minLInd ]) [5 S4 A32])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 453 452 454 61 ("./CppStatUtilities.cc") 520)

(insn 454 453 455 61 ./CppStatUtilities.cc:520 (set (reg/v:DF 102 [ minVal ])
        (mem:DF (reg/v/f:DI 116 [ x ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 455 454 456 61 ./CppStatUtilities.cc:520 (set (reg/v:DF 101 [ maxVal ])
        (reg/v:DF 102 [ minVal ])) -1 (nil)
    (nil))

(insn 456 455 457 61 ./CppStatUtilities.cc:520 (set (reg/v:SI 103 [ count ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 457 456 458 61 ./CppStatUtilities.cc:520 (set (reg/v:SI 84 [ i.502 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 458 457 459 61 ./CppStatUtilities.cc:520 (set (reg:DI 98 [ prephitmp.418 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 459 458 460 61 ./CppStatUtilities.cc:520 (set (reg:QI 97 [ prephitmp.420 ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(insn 460 459 461 61 ./CppStatUtilities.cc:520 (set (reg:SI 100 [ prephitmp.406 ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(insn 461 460 462 61 ./CppStatUtilities.cc:520 (set (reg:DI 99 [ prephitmp.410 ])
        (const_int 48 [0x30])) -1 (nil)
    (nil))

(jump_insn 462 461 463 61 ./CppStatUtilities.cc:520 (set (pc)
        (label_ref 339)) -1 (nil)
    (nil))
;; End of basic block 61, registers live:
 (nil)

(barrier 463 462 464)

;; Start of basic block 62, registers live: (nil)
(code_label 464 463 465 62 169 "" [1 uses])

(note 465 464 466 62 [bb 62] NOTE_INSN_BASIC_BLOCK)

(note 466 465 467 62 ("./CppStatUtilities.cc") 519)

(insn 467 466 468 62 ./CppStatUtilities.cc:519 (set (mem:SI (reg/v/f:DI 107 [ minLInd ]) [5 S4 A32])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(note 468 467 469 62 ("./CppStatUtilities.cc") 520)

(insn 469 468 470 62 ./CppStatUtilities.cc:520 (set (reg/v:DF 102 [ minVal ])
        (mem:DF (reg/v/f:DI 116 [ x ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 470 469 471 62 ./CppStatUtilities.cc:520 (set (reg/v:DF 101 [ maxVal ])
        (reg/v:DF 102 [ minVal ])) -1 (nil)
    (nil))

(insn 471 470 472 62 ./CppStatUtilities.cc:520 (set (reg/v:SI 85 [ i.501 ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))
;; End of basic block 62, registers live:
 (nil)

;; Start of basic block 63, registers live: (nil)
(code_label 472 471 473 63 180 "" [1 uses])

(note 473 472 474 63 [bb 63] NOTE_INSN_BASIC_BLOCK)

(insn 474 473 475 63 ./CppStatUtilities.cc:520 (set (reg/v:SI 103 [ count ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 475 474 476 63 ./CppStatUtilities.cc:520 (set (pc)
        (label_ref 162)) -1 (nil)
    (nil))
;; End of basic block 63, registers live:
 (nil)

(barrier 476 475 477)

;; Start of basic block 64, registers live: (nil)
(code_label 477 476 478 64 178 "" [1 uses])

(note 478 477 479 64 [bb 64] NOTE_INSN_BASIC_BLOCK)

(note 479 478 480 64 ("./CppStatUtilities.cc") 570)

(insn 480 479 481 64 ./CppStatUtilities.cc:570 (parallel [
            (set (reg/v:SI 84 [ i.502 ])
                (plus:SI (reg/v:SI 85 [ i.501 ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 481 480 482 64 ./CppStatUtilities.cc:570 (set (reg/v:DF 101 [ maxVal ])
        (reg/v:DF 102 [ minVal ])) -1 (nil)
    (nil))

(insn 482 481 483 64 ./CppStatUtilities.cc:570 (set (reg:DI 98 [ prephitmp.418 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 483 482 484 64 ./CppStatUtilities.cc:570 (set (reg:SI 100 [ prephitmp.406 ])
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(insn 484 483 485 64 ./CppStatUtilities.cc:570 (set (reg:DI 99 [ prephitmp.410 ])
        (const_int 48 [0x30])) -1 (nil)
    (nil))

(jump_insn 485 484 486 64 ./CppStatUtilities.cc:570 (set (pc)
        (label_ref 356)) -1 (nil)
    (nil))
;; End of basic block 64, registers live:
 (nil)

(barrier 486 485 487)

(note 487 486 488 NOTE_INSN_FUNCTION_END)

(note 488 487 0 ("./CppStatUtilities.cc") 601)


;; Function void hist(double*, int, double*, int, int, double**, double**, double**) (_Z4histPdiS_iiPS_S0_S0_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Edge 0->3 redirected to 4
Forwarding edge 0->2 to 5 failed.
Forwarding edge 0->2 to 5 failed.
Deleting block 3.
Edge 6->8 redirected to 9
Forwarding edge 6->7 to 10 failed.
Forwarding edge 6->7 to 10 failed.
Deleting block 8.
Deleted label in block 12.
Deleted label in block 15.
Deleted label in block 18.
Deleted label in block 21.
Merged 23 and 24 without moving.
Merged 23 and 25 without moving.


try_optimize_cfg iteration 2

Forwarding edge 0->2 to 5 failed.
Forwarding edge 6->7 to 10 failed.


try_optimize_cfg iteration 1

Forwarding edge 0->1 to 3 failed.
Forwarding edge 4->5 to 7 failed.
(note 1 0 13 ("./CppStatUtilities.cc") 449)

;; Start of basic block 0, registers live: (nil)
(note 13 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 13 4 0 ./CppStatUtilities.cc:449 (set (reg/v/f:DI 87 [ x ])
        (reg:DI 5 di [ x ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:449 (set (reg/v:SI 88 [ xSize ])
        (reg:SI 4 si [ xSize ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:449 (set (reg/v/f:DI 89 [ y ])
        (reg:DI 1 dx [ y ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:449 (set (reg/v:SI 90 [ ySize ])
        (reg:SI 2 cx [ ySize ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:449 (set (reg/v:SI 91 [ nBins ])
        (reg:SI 37 r8 [ nBins ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppStatUtilities.cc:449 (set (reg/v/f:DI 92 [ _hmid ])
        (reg:DI 38 r9 [ _hmid ])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppStatUtilities.cc:449 (set (reg/v/f:DI 93 [ _hx ])
        (mem/f/c/i:DI (reg/f:DI 53 virtual-incoming-args) [23 _hx+0 S8 A64])) -1 (nil)
    (expr_list:REG_EQUIV (mem/f/c/i:DI (reg/f:DI 53 virtual-incoming-args) [23 _hx+0 S8 A64])
        (nil)))

(insn 10 9 11 0 ./CppStatUtilities.cc:449 (set (reg/v/f:DI 94 [ _hy ])
        (mem/f/c/i:DI (plus:DI (reg/f:DI 53 virtual-incoming-args)
                (const_int 8 [0x8])) [23 _hy+0 S8 A64])) -1 (nil)
    (expr_list:REG_EQUIV (mem/f/c/i:DI (plus:DI (reg/f:DI 53 virtual-incoming-args)
                (const_int 8 [0x8])) [23 _hy+0 S8 A64])
        (nil)))

(note 11 10 15 0 NOTE_INSN_FUNCTION_BEG)

(note 15 11 16 0 ("./CppStatUtilities.cc") 451)

(insn 16 15 17 0 ./CppStatUtilities.cc:451 (set (reg:SI 4 si)
        (reg/v:SI 88 [ xSize ])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:451 (set (reg:DI 5 di)
        (reg/v/f:DI 87 [ x ])) -1 (nil)
    (nil))

(call_insn 18 17 19 0 ./CppStatUtilities.cc:451 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("min") [flags 0x41] <function_decl 0x2b5eb631ee00 min>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 19 18 20 0 ./CppStatUtilities.cc:451 (set (reg/v:DF 79 [ xmin ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 20 19 21 0 ("./CppStatUtilities.cc") 452)

(insn 21 20 22 0 ./CppStatUtilities.cc:452 (set (reg:SI 4 si)
        (reg/v:SI 88 [ xSize ])) -1 (nil)
    (nil))

(insn 22 21 23 0 ./CppStatUtilities.cc:452 (set (reg:DI 5 di)
        (reg/v/f:DI 87 [ x ])) -1 (nil)
    (nil))

(call_insn 23 22 24 0 ./CppStatUtilities.cc:452 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("max") [flags 0x41] <function_decl 0x2b5eb631ef00 max>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 24 23 25 0 ./CppStatUtilities.cc:452 (set (reg/v:DF 78 [ xmax ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 25 24 26 0 ("./CppStatUtilities.cc") 454)

(insn 26 25 27 0 ./CppStatUtilities.cc:454 (set (reg:SI 4 si)
        (reg/v:SI 90 [ ySize ])) -1 (nil)
    (nil))

(insn 27 26 28 0 ./CppStatUtilities.cc:454 (set (reg:DI 5 di)
        (reg/v/f:DI 89 [ y ])) -1 (nil)
    (nil))

(call_insn 28 27 29 0 ./CppStatUtilities.cc:454 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("min") [flags 0x41] <function_decl 0x2b5eb631ee00 min>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 29 28 30 0 ./CppStatUtilities.cc:454 (set (reg/v:DF 77 [ ymin ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 30 29 31 0 ("./CppStatUtilities.cc") 455)

(insn 31 30 32 0 ./CppStatUtilities.cc:455 (set (reg:SI 4 si)
        (reg/v:SI 90 [ ySize ])) -1 (nil)
    (nil))

(insn 32 31 33 0 ./CppStatUtilities.cc:455 (set (reg:DI 5 di)
        (reg/v/f:DI 89 [ y ])) -1 (nil)
    (nil))

(call_insn 33 32 34 0 ./CppStatUtilities.cc:455 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("max") [flags 0x41] <function_decl 0x2b5eb631ef00 max>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 34 33 35 0 ./CppStatUtilities.cc:455 (set (reg/v:DF 76 [ ymax ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 35 34 36 0 ("./CppStatUtilities.cc") 457)

(insn 36 35 37 0 ./CppStatUtilities.cc:457 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 77 [ ymin ])
            (reg/v:DF 79 [ xmin ]))) -1 (nil)
    (nil))

(jump_insn 37 36 223 0 ./CppStatUtilities.cc:457 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 41)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 223 37 38 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(jump_insn 38 223 39 1 ./CppStatUtilities.cc:457 (set (pc)
        (label_ref 46)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 39 38 41)

;; Start of basic block 2, registers live: (nil)
(code_label 41 39 42 2 331 "" [1 uses])

(note 42 41 43 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 43 42 44 2 ./CppStatUtilities.cc:457 (set (reg/v:DF 75 [ tmin ])
        (reg/v:DF 79 [ xmin ])) -1 (nil)
    (nil))

(jump_insn 44 43 45 2 ./CppStatUtilities.cc:457 (set (pc)
        (label_ref 49)) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

(barrier 45 44 46)

;; Start of basic block 3, registers live: (nil)
(code_label 46 45 47 3 329 "" [1 uses])

(note 47 46 48 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 48 47 49 3 ./CppStatUtilities.cc:457 (set (reg/v:DF 75 [ tmin ])
        (reg/v:DF 77 [ ymin ])) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 49 48 50 4 332 "" [1 uses])

(note 50 49 51 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 51 50 52 4 ("./CppStatUtilities.cc") 458)

(insn 52 51 53 4 ./CppStatUtilities.cc:458 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg/v:DF 78 [ xmax ])
            (reg/v:DF 76 [ ymax ]))) -1 (nil)
    (nil))

(jump_insn 53 52 225 4 ./CppStatUtilities.cc:458 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref:DI 57)
            (pc))) 533 {*jcc_1} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 225 53 54 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(jump_insn 54 225 55 5 ./CppStatUtilities.cc:458 (set (pc)
        (label_ref 62)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 55 54 57)

;; Start of basic block 6, registers live: (nil)
(code_label 57 55 58 6 335 "" [1 uses])

(note 58 57 59 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 59 58 60 6 ./CppStatUtilities.cc:458 (set (reg/v:DF 74 [ tmax ])
        (reg/v:DF 78 [ xmax ])) -1 (nil)
    (nil))

(jump_insn 60 59 61 6 ./CppStatUtilities.cc:458 (set (pc)
        (label_ref 65)) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

(barrier 61 60 62)

;; Start of basic block 7, registers live: (nil)
(code_label 62 61 63 7 333 "" [1 uses])

(note 63 62 64 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 64 63 65 7 ./CppStatUtilities.cc:458 (set (reg/v:DF 74 [ tmax ])
        (reg/v:DF 76 [ ymax ])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 65 64 66 8 336 "" [1 uses])

(note 66 65 67 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 67 66 68 8 ("./CppStatUtilities.cc") 460)

(insn 68 67 69 8 ./CppStatUtilities.cc:460 (set (reg:DI 86 [ D.35109 ])
        (sign_extend:DI (reg/v:SI 91 [ nBins ]))) -1 (nil)
    (nil))

(insn 69 68 70 8 ./CppStatUtilities.cc:460 (parallel [
            (set (reg:DI 95)
                (ashift:DI (reg:DI 86 [ D.35109 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 70 69 71 8 ./CppStatUtilities.cc:460 (set (reg:DI 5 di)
        (reg:DI 95)) -1 (nil)
    (nil))

(call_insn 71 70 72 8 ./CppStatUtilities.cc:460 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 72 71 73 8 ./CppStatUtilities.cc:460 (set (reg/f:DI 96)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 96)
        (nil)))

(insn 73 72 74 8 ./CppStatUtilities.cc:460 (set (reg:DI 58 [ ivtmp.606 ])
        (reg/f:DI 96)) -1 (nil)
    (nil))

(insn 74 73 75 8 ./CppStatUtilities.cc:460 (set (reg/v/f:DI 73 [ hmid ])
        (reg:DI 58 [ ivtmp.606 ])) -1 (nil)
    (nil))

(note 75 74 76 8 ("./CppStatUtilities.cc") 461)

(insn 76 75 77 8 ./CppStatUtilities.cc:461 (set (reg:DI 4 si)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(insn 77 76 78 8 ./CppStatUtilities.cc:461 (set (reg:DI 5 di)
        (reg:DI 86 [ D.35109 ])) -1 (nil)
    (nil))

(call_insn 78 77 79 8 ./CppStatUtilities.cc:461 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b5eb3cbd200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 79 78 80 8 ./CppStatUtilities.cc:461 (set (reg/f:DI 97)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 97)
        (nil)))

(insn 80 79 81 8 ./CppStatUtilities.cc:461 (set (reg:DI 85 [ D.35113 ])
        (reg/f:DI 97)) -1 (nil)
    (nil))

(insn 81 80 82 8 ./CppStatUtilities.cc:461 (set (reg/v/f:DI 72 [ hx ])
        (reg:DI 85 [ D.35113 ])) -1 (nil)
    (nil))

(note 82 81 83 8 ("./CppStatUtilities.cc") 462)

(insn 83 82 84 8 ./CppStatUtilities.cc:462 (set (reg:DI 4 si)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(insn 84 83 85 8 ./CppStatUtilities.cc:462 (set (reg:DI 5 di)
        (reg:DI 86 [ D.35109 ])) -1 (nil)
    (nil))

(call_insn 85 84 86 8 ./CppStatUtilities.cc:462 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b5eb3cbd200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 86 85 87 8 ./CppStatUtilities.cc:462 (set (reg/f:DI 98)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 98)
        (nil)))

(insn 87 86 88 8 ./CppStatUtilities.cc:462 (set (reg:DI 84 [ D.35114 ])
        (reg/f:DI 98)) -1 (nil)
    (nil))

(insn 88 87 89 8 ./CppStatUtilities.cc:462 (set (reg/v/f:DI 71 [ hy ])
        (reg:DI 84 [ D.35114 ])) -1 (nil)
    (nil))

(note 89 88 90 8 ("./CppStatUtilities.cc") 464)

(insn 90 89 91 8 ./CppStatUtilities.cc:464 (set (reg:DF 99)
        (minus:DF (reg/v:DF 74 [ tmax ])
            (reg/v:DF 75 [ tmin ]))) -1 (nil)
    (nil))

(insn 91 90 92 8 ./CppStatUtilities.cc:464 (set (reg:DF 100)
        (float:DF (reg/v:SI 91 [ nBins ]))) -1 (nil)
    (nil))

(insn 92 91 93 8 ./CppStatUtilities.cc:464 (set (reg/v:DF 70 [ dt ])
        (div:DF (reg:DF 99)
            (reg:DF 100))) -1 (nil)
    (nil))

(note 93 92 94 8 ("./CppStatUtilities.cc") 469)

(insn 94 93 95 8 ./CppStatUtilities.cc:469 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 91 [ nBins ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 95 94 97 8 ./CppStatUtilities.cc:469 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 114)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(note 97 95 98 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 98 97 99 9 ./CppStatUtilities.cc:469 (set (reg:DF 101)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC5") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 99 98 100 9 ./CppStatUtilities.cc:469 (set (reg:DF 102)
        (mult:DF (reg/v:DF 70 [ dt ])
            (reg:DF 101))) -1 (nil)
    (nil))

(insn 100 99 101 9 ./CppStatUtilities.cc:469 (set (reg:DF 63 [ pretmp.572 ])
        (plus:DF (reg/v:DF 75 [ tmin ])
            (reg:DF 102))) -1 (nil)
    (nil))

(insn 101 100 102 9 ./CppStatUtilities.cc:469 (set (reg/v:SI 67 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 102 101 103 10 339 "" [1 uses])

(note 103 102 104 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 104 103 105 10 ("./CppStatUtilities.cc") 470)

(insn 105 104 106 10 ./CppStatUtilities.cc:470 (set (reg:DF 103)
        (float:DF (reg/v:SI 67 [ i ]))) -1 (nil)
    (nil))

(insn 106 105 107 10 ./CppStatUtilities.cc:470 (set (reg:DF 104)
        (mult:DF (reg:DF 103)
            (reg/v:DF 70 [ dt ]))) -1 (nil)
    (nil))

(insn 107 106 108 10 ./CppStatUtilities.cc:470 (set (reg:DF 105)
        (plus:DF (reg:DF 63 [ pretmp.572 ])
            (reg:DF 104))) -1 (nil)
    (nil))

(insn 108 107 109 10 ./CppStatUtilities.cc:470 (set (mem:DF (reg:DI 58 [ ivtmp.606 ]) [3 S8 A64])
        (reg:DF 105)) -1 (nil)
    (nil))

(note 109 108 110 10 ("./CppStatUtilities.cc") 469)

(insn 110 109 111 10 ./CppStatUtilities.cc:469 (parallel [
            (set (reg/v:SI 67 [ i ])
                (plus:SI (reg/v:SI 67 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 111 110 112 10 ./CppStatUtilities.cc:469 (parallel [
            (set (reg:DI 58 [ ivtmp.606 ])
                (plus:DI (reg:DI 58 [ ivtmp.606 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 112 111 113 10 ./CppStatUtilities.cc:469 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 67 [ i ])
            (reg/v:SI 91 [ nBins ]))) -1 (nil)
    (nil))

(jump_insn 113 112 114 10 ./CppStatUtilities.cc:469 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 102)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 10, registers live:
 (nil)

;; Start of basic block 11, registers live: (nil)
(code_label 114 113 115 11 337 "" [1 uses])

(note 115 114 116 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(note 116 115 117 11 ("./CppStatUtilities.cc") 473)

(insn 117 116 118 11 ./CppStatUtilities.cc:473 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 88 [ xSize ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 118 117 120 11 ./CppStatUtilities.cc:473 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 145)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(note 120 118 121 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(insn 121 120 122 12 ./CppStatUtilities.cc:473 (set (reg:DI 59 [ ivtmp.601 ])
        (reg/v/f:DI 87 [ x ])) -1 (nil)
    (nil))

(insn 122 121 123 12 ./CppStatUtilities.cc:473 (set (reg/v:SI 66 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(code_label 123 122 124 13 342 "" [1 uses])

(note 124 123 125 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(note 125 124 126 13 ("./CppStatUtilities.cc") 474)

(insn 126 125 127 13 ./CppStatUtilities.cc:474 (set (reg:DF 107)
        (mem:DF (reg:DI 59 [ ivtmp.601 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 127 126 128 13 ./CppStatUtilities.cc:474 (set (reg:DF 106)
        (minus:DF (reg:DF 107)
            (reg/v:DF 75 [ tmin ]))) -1 (nil)
    (nil))

(insn 128 127 129 13 ./CppStatUtilities.cc:474 (set (reg:DF 108)
        (div:DF (reg:DF 106)
            (reg/v:DF 70 [ dt ]))) -1 (nil)
    (nil))

(insn 129 128 130 13 ./CppStatUtilities.cc:474 (set (reg:DF 21 xmm0)
        (reg:DF 108)) -1 (nil)
    (nil))

(call_insn/u 130 129 131 13 ./CppStatUtilities.cc:474 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("floor") [flags 0x41] <function_decl 0x2b5eb3c39300 floor>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 131 130 132 13 ./CppStatUtilities.cc:474 (set (reg:DF 83 [ D.35140 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 132 131 133 13 ./CppStatUtilities.cc:474 (set (reg:SI 109)
        (fix:SI (reg:DF 83 [ D.35140 ]))) -1 (nil)
    (nil))

(insn 133 132 134 13 ./CppStatUtilities.cc:474 (set (reg:DI 110)
        (sign_extend:DI (reg:SI 109))) -1 (nil)
    (nil))

(insn 134 133 135 13 ./CppStatUtilities.cc:474 (parallel [
            (set (reg:DI 111)
                (ashift:DI (reg:DI 110)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 135 134 136 13 ./CppStatUtilities.cc:474 (parallel [
            (set (reg:DI 82 [ D.35145 ])
                (plus:DI (reg:DI 111)
                    (reg/v/f:DI 72 [ hx ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 136 135 137 13 ./CppStatUtilities.cc:474 (set (reg:DF 112)
        (mem:DF (reg:DI 82 [ D.35145 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 137 136 138 13 ./CppStatUtilities.cc:474 (set (reg:DF 113)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC6") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 138 137 139 13 ./CppStatUtilities.cc:474 (set (reg:DF 114)
        (plus:DF (reg:DF 112)
            (reg:DF 113))) -1 (nil)
    (nil))

(insn 139 138 140 13 ./CppStatUtilities.cc:474 (set (mem:DF (reg:DI 82 [ D.35145 ]) [3 S8 A64])
        (reg:DF 114)) -1 (nil)
    (nil))

(note 140 139 141 13 ("./CppStatUtilities.cc") 473)

(insn 141 140 142 13 ./CppStatUtilities.cc:473 (parallel [
            (set (reg/v:SI 66 [ i ])
                (plus:SI (reg/v:SI 66 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 142 141 143 13 ./CppStatUtilities.cc:473 (parallel [
            (set (reg:DI 59 [ ivtmp.601 ])
                (plus:DI (reg:DI 59 [ ivtmp.601 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 143 142 144 13 ./CppStatUtilities.cc:473 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 66 [ i ])
            (reg/v:SI 88 [ xSize ]))) -1 (nil)
    (nil))

(jump_insn 144 143 145 13 ./CppStatUtilities.cc:473 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 123)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 13, registers live:
 (nil)

;; Start of basic block 14, registers live: (nil)
(code_label 145 144 146 14 340 "" [1 uses])

(note 146 145 147 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(note 147 146 148 14 ("./CppStatUtilities.cc") 477)

(insn 148 147 149 14 ./CppStatUtilities.cc:477 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 90 [ ySize ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 149 148 151 14 ./CppStatUtilities.cc:477 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 176)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 14, registers live:
 (nil)

;; Start of basic block 15, registers live: (nil)
(note 151 149 152 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(insn 152 151 153 15 ./CppStatUtilities.cc:477 (set (reg:DI 60 [ ivtmp.596 ])
        (reg/v/f:DI 89 [ y ])) -1 (nil)
    (nil))

(insn 153 152 154 15 ./CppStatUtilities.cc:477 (set (reg/v:SI 65 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(code_label 154 153 155 16 345 "" [1 uses])

(note 155 154 156 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(note 156 155 157 16 ("./CppStatUtilities.cc") 478)

(insn 157 156 158 16 ./CppStatUtilities.cc:478 (set (reg:DF 116)
        (mem:DF (reg:DI 60 [ ivtmp.596 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 158 157 159 16 ./CppStatUtilities.cc:478 (set (reg:DF 115)
        (minus:DF (reg:DF 116)
            (reg/v:DF 75 [ tmin ]))) -1 (nil)
    (nil))

(insn 159 158 160 16 ./CppStatUtilities.cc:478 (set (reg:DF 117)
        (div:DF (reg:DF 115)
            (reg/v:DF 70 [ dt ]))) -1 (nil)
    (nil))

(insn 160 159 161 16 ./CppStatUtilities.cc:478 (set (reg:DF 21 xmm0)
        (reg:DF 117)) -1 (nil)
    (nil))

(call_insn/u 161 160 162 16 ./CppStatUtilities.cc:478 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("floor") [flags 0x41] <function_decl 0x2b5eb3c39300 floor>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 162 161 163 16 ./CppStatUtilities.cc:478 (set (reg:DF 81 [ D.35159 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 163 162 164 16 ./CppStatUtilities.cc:478 (set (reg:SI 118)
        (fix:SI (reg:DF 81 [ D.35159 ]))) -1 (nil)
    (nil))

(insn 164 163 165 16 ./CppStatUtilities.cc:478 (set (reg:DI 119)
        (sign_extend:DI (reg:SI 118))) -1 (nil)
    (nil))

(insn 165 164 166 16 ./CppStatUtilities.cc:478 (parallel [
            (set (reg:DI 120)
                (ashift:DI (reg:DI 119)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 166 165 167 16 ./CppStatUtilities.cc:478 (parallel [
            (set (reg:DI 80 [ D.35164 ])
                (plus:DI (reg:DI 120)
                    (reg/v/f:DI 71 [ hy ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 167 166 168 16 ./CppStatUtilities.cc:478 (set (reg:DF 121)
        (mem:DF (reg:DI 80 [ D.35164 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 168 167 169 16 ./CppStatUtilities.cc:478 (set (reg:DF 122)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC6") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 169 168 170 16 ./CppStatUtilities.cc:478 (set (reg:DF 123)
        (plus:DF (reg:DF 121)
            (reg:DF 122))) -1 (nil)
    (nil))

(insn 170 169 171 16 ./CppStatUtilities.cc:478 (set (mem:DF (reg:DI 80 [ D.35164 ]) [3 S8 A64])
        (reg:DF 123)) -1 (nil)
    (nil))

(note 171 170 172 16 ("./CppStatUtilities.cc") 477)

(insn 172 171 173 16 ./CppStatUtilities.cc:477 (parallel [
            (set (reg/v:SI 65 [ i ])
                (plus:SI (reg/v:SI 65 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 173 172 174 16 ./CppStatUtilities.cc:477 (parallel [
            (set (reg:DI 60 [ ivtmp.596 ])
                (plus:DI (reg:DI 60 [ ivtmp.596 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 174 173 175 16 ./CppStatUtilities.cc:477 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 65 [ i ])
            (reg/v:SI 90 [ ySize ]))) -1 (nil)
    (nil))

(jump_insn 175 174 176 16 ./CppStatUtilities.cc:477 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 154)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(code_label 176 175 177 17 343 "" [1 uses])

(note 177 176 178 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(note 178 177 179 17 ("./CppStatUtilities.cc") 482)

(insn 179 178 180 17 ./CppStatUtilities.cc:482 (set (reg:DF 124)
        (float:DF (reg/v:SI 88 [ xSize ]))) -1 (nil)
    (nil))

(insn 180 179 181 17 ./CppStatUtilities.cc:482 (set (reg/v:DF 69 [ fx ])
        (mult:DF (reg:DF 124)
            (reg/v:DF 70 [ dt ]))) -1 (nil)
    (nil))

(note 181 180 182 17 ("./CppStatUtilities.cc") 483)

(insn 182 181 183 17 ./CppStatUtilities.cc:483 (set (reg:DF 125)
        (float:DF (reg/v:SI 90 [ ySize ]))) -1 (nil)
    (nil))

(insn 183 182 184 17 ./CppStatUtilities.cc:483 (set (reg/v:DF 68 [ fy ])
        (mult:DF (reg:DF 125)
            (reg/v:DF 70 [ dt ]))) -1 (nil)
    (nil))

(note 184 183 185 17 ("./CppStatUtilities.cc") 485)

(insn 185 184 186 17 ./CppStatUtilities.cc:485 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 91 [ nBins ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 186 185 188 17 ./CppStatUtilities.cc:485 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 208)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(note 188 186 189 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(insn 189 188 190 18 ./CppStatUtilities.cc:485 (set (reg/v:SI 64 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 190 189 191 18 ./CppStatUtilities.cc:485 (set (reg:DI 62 [ ivtmp.588 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

;; Start of basic block 19, registers live: (nil)
(code_label 191 190 192 19 348 "" [1 uses])

(note 192 191 193 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(note 193 192 194 19 ("./CppStatUtilities.cc") 449)

(insn 194 193 195 19 ./CppStatUtilities.cc:449 (set (reg:DI 61 [ D.36425 ])
        (reg:DI 62 [ ivtmp.588 ])) -1 (nil)
    (nil))

(note 195 194 196 19 ("./CppStatUtilities.cc") 487)

(insn 196 195 197 19 ./CppStatUtilities.cc:487 (set (reg:DF 126)
        (mem:DF (plus:DI (mult:DI (reg:DI 61 [ D.36425 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 72 [ hx ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 197 196 198 19 ./CppStatUtilities.cc:487 (set (reg:DF 127)
        (div:DF (reg:DF 126)
            (reg/v:DF 69 [ fx ]))) -1 (nil)
    (nil))

(insn 198 197 199 19 ./CppStatUtilities.cc:487 (set (mem:DF (plus:DI (mult:DI (reg:DI 61 [ D.36425 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 72 [ hx ])) [3 S8 A64])
        (reg:DF 127)) -1 (nil)
    (nil))

(note 199 198 200 19 ("./CppStatUtilities.cc") 488)

(insn 200 199 201 19 ./CppStatUtilities.cc:488 (set (reg:DF 128)
        (mem:DF (plus:DI (mult:DI (reg:DI 61 [ D.36425 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 71 [ hy ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 201 200 202 19 ./CppStatUtilities.cc:488 (set (reg:DF 129)
        (div:DF (reg:DF 128)
            (reg/v:DF 68 [ fy ]))) -1 (nil)
    (nil))

(insn 202 201 203 19 ./CppStatUtilities.cc:488 (set (mem:DF (plus:DI (mult:DI (reg:DI 61 [ D.36425 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 71 [ hy ])) [3 S8 A64])
        (reg:DF 129)) -1 (nil)
    (nil))

(note 203 202 204 19 ("./CppStatUtilities.cc") 485)

(insn 204 203 205 19 ./CppStatUtilities.cc:485 (parallel [
            (set (reg/v:SI 64 [ i ])
                (plus:SI (reg/v:SI 64 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 205 204 206 19 ./CppStatUtilities.cc:485 (parallel [
            (set (reg:DI 62 [ ivtmp.588 ])
                (plus:DI (reg:DI 62 [ ivtmp.588 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 206 205 207 19 ./CppStatUtilities.cc:485 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 64 [ i ])
            (reg/v:SI 91 [ nBins ]))) -1 (nil)
    (nil))

(jump_insn 207 206 208 19 ./CppStatUtilities.cc:485 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 191)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 19, registers live:
 (nil)

;; Start of basic block 20, registers live: (nil)
(code_label 208 207 209 20 346 "" [1 uses])

(note 209 208 210 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(note 210 209 211 20 ("./CppStatUtilities.cc") 491)

(insn 211 210 212 20 ./CppStatUtilities.cc:491 (set (mem/f:DI (reg/v/f:DI 92 [ _hmid ]) [22 S8 A64])
        (reg/v/f:DI 73 [ hmid ])) -1 (nil)
    (nil))

(note 212 211 213 20 ("./CppStatUtilities.cc") 492)

(insn 213 212 214 20 ./CppStatUtilities.cc:492 (set (mem/f:DI (reg/v/f:DI 93 [ _hx ]) [22 S8 A64])
        (reg/v/f:DI 72 [ hx ])) -1 (nil)
    (nil))

(note 214 213 215 20 ("./CppStatUtilities.cc") 493)

(insn 215 214 216 20 ./CppStatUtilities.cc:493 (set (mem/f:DI (reg/v/f:DI 94 [ _hy ]) [22 S8 A64])
        (reg/v/f:DI 71 [ hy ])) -1 (nil)
    (nil))
;; End of basic block 20, registers live:
 (nil)

(note 216 215 217 NOTE_INSN_FUNCTION_END)

(note 217 216 0 ("./CppStatUtilities.cc") 494)


;; Function void hist(double*, int, int, double**, double**) (_Z4histPdiiPS_S0_)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 5.
Deleted label in block 8.
Merged 10 and 11 without moving.
Merged 10 and 12 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 10 ("./CppStatUtilities.cc") 412)

;; Start of basic block 0, registers live: (nil)
(note 10 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 10 4 0 ./CppStatUtilities.cc:412 (set (reg/v/f:DI 75 [ x ])
        (reg:DI 5 di [ x ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:412 (set (reg/v:SI 76 [ xSize ])
        (reg:SI 4 si [ xSize ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:412 (set (reg/v:SI 77 [ nBins ])
        (reg:SI 1 dx [ nBins ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:412 (set (reg/v/f:DI 78 [ _hx ])
        (reg:DI 2 cx [ _hx ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:412 (set (reg/v/f:DI 79 [ _hy ])
        (reg:DI 37 r8 [ _hy ])) -1 (nil)
    (nil))

(note 8 7 12 0 NOTE_INSN_FUNCTION_BEG)

(note 12 8 13 0 ("./CppStatUtilities.cc") 414)

(insn 13 12 14 0 ./CppStatUtilities.cc:414 (set (reg:SI 4 si)
        (reg/v:SI 76 [ xSize ])) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppStatUtilities.cc:414 (set (reg:DI 5 di)
        (reg/v/f:DI 75 [ x ])) -1 (nil)
    (nil))

(call_insn 15 14 16 0 ./CppStatUtilities.cc:414 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("min") [flags 0x41] <function_decl 0x2b5eb631ee00 min>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 16 15 17 0 ./CppStatUtilities.cc:414 (set (reg/v:DF 71 [ xmin ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 17 16 18 0 ("./CppStatUtilities.cc") 415)

(insn 18 17 19 0 ./CppStatUtilities.cc:415 (set (reg:SI 4 si)
        (reg/v:SI 76 [ xSize ])) -1 (nil)
    (nil))

(insn 19 18 20 0 ./CppStatUtilities.cc:415 (set (reg:DI 5 di)
        (reg/v/f:DI 75 [ x ])) -1 (nil)
    (nil))

(call_insn 20 19 21 0 ./CppStatUtilities.cc:415 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("max") [flags 0x41] <function_decl 0x2b5eb631ef00 max>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 21 20 22 0 ./CppStatUtilities.cc:415 (set (reg/v:DF 70 [ xmax ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 22 21 23 0 ("./CppStatUtilities.cc") 417)

(insn 23 22 24 0 ./CppStatUtilities.cc:417 (set (reg:DI 74 [ D.35027 ])
        (sign_extend:DI (reg/v:SI 77 [ nBins ]))) -1 (nil)
    (nil))

(insn 24 23 25 0 ./CppStatUtilities.cc:417 (parallel [
            (set (reg:DI 80)
                (ashift:DI (reg:DI 74 [ D.35027 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 25 24 26 0 ./CppStatUtilities.cc:417 (set (reg:DI 5 di)
        (reg:DI 80)) -1 (nil)
    (nil))

(call_insn 26 25 27 0 ./CppStatUtilities.cc:417 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 27 26 28 0 ./CppStatUtilities.cc:417 (set (reg/f:DI 81)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 81)
        (nil)))

(insn 28 27 29 0 ./CppStatUtilities.cc:417 (set (reg:DI 58 [ ivtmp.664 ])
        (reg/f:DI 81)) -1 (nil)
    (nil))

(insn 29 28 30 0 ./CppStatUtilities.cc:417 (set (reg/v/f:DI 69 [ hx ])
        (reg:DI 58 [ ivtmp.664 ])) -1 (nil)
    (nil))

(note 30 29 31 0 ("./CppStatUtilities.cc") 418)

(insn 31 30 32 0 ./CppStatUtilities.cc:418 (set (reg:DI 4 si)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(insn 32 31 33 0 ./CppStatUtilities.cc:418 (set (reg:DI 5 di)
        (reg:DI 74 [ D.35027 ])) -1 (nil)
    (nil))

(call_insn 33 32 34 0 ./CppStatUtilities.cc:418 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("calloc") [flags 0x41] <function_decl 0x2b5eb3cbd200 calloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 34 33 35 0 ./CppStatUtilities.cc:418 (set (reg/f:DI 82)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 82)
        (nil)))

(insn 35 34 36 0 ./CppStatUtilities.cc:418 (set (reg:DI 61 [ ivtmp.654 ])
        (reg/f:DI 82)) -1 (nil)
    (nil))

(insn 36 35 37 0 ./CppStatUtilities.cc:418 (set (reg/v/f:DI 68 [ hy ])
        (reg:DI 61 [ ivtmp.654 ])) -1 (nil)
    (nil))

(note 37 36 38 0 ("./CppStatUtilities.cc") 420)

(insn 38 37 39 0 ./CppStatUtilities.cc:420 (set (reg:DF 83)
        (minus:DF (reg/v:DF 70 [ xmax ])
            (reg/v:DF 71 [ xmin ]))) -1 (nil)
    (nil))

(insn 39 38 40 0 ./CppStatUtilities.cc:420 (set (reg:DF 84)
        (float:DF (reg/v:SI 77 [ nBins ]))) -1 (nil)
    (nil))

(insn 40 39 41 0 ./CppStatUtilities.cc:420 (set (reg/v:DF 67 [ dx ])
        (div:DF (reg:DF 83)
            (reg:DF 84))) -1 (nil)
    (nil))

(note 41 40 42 0 ("./CppStatUtilities.cc") 425)

(insn 42 41 43 0 ./CppStatUtilities.cc:425 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 77 [ nBins ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 43 42 45 0 ./CppStatUtilities.cc:425 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 62)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 45 43 46 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 46 45 47 1 ./CppStatUtilities.cc:425 (set (reg:DF 85)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC7") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 47 46 48 1 ./CppStatUtilities.cc:425 (set (reg:DF 86)
        (mult:DF (reg/v:DF 67 [ dx ])
            (reg:DF 85))) -1 (nil)
    (nil))

(insn 48 47 49 1 ./CppStatUtilities.cc:425 (set (reg:DF 62 [ pretmp.640 ])
        (plus:DF (reg/v:DF 71 [ xmin ])
            (reg:DF 86))) -1 (nil)
    (nil))

(insn 49 48 50 1 ./CppStatUtilities.cc:425 (set (reg/v:SI 65 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 50 49 51 2 453 "" [1 uses])

(note 51 50 52 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 52 51 53 2 ("./CppStatUtilities.cc") 426)

(insn 53 52 54 2 ./CppStatUtilities.cc:426 (set (reg:DF 87)
        (float:DF (reg/v:SI 65 [ i ]))) -1 (nil)
    (nil))

(insn 54 53 55 2 ./CppStatUtilities.cc:426 (set (reg:DF 88)
        (mult:DF (reg:DF 87)
            (reg/v:DF 67 [ dx ]))) -1 (nil)
    (nil))

(insn 55 54 56 2 ./CppStatUtilities.cc:426 (set (reg:DF 89)
        (plus:DF (reg:DF 88)
            (reg:DF 62 [ pretmp.640 ]))) -1 (nil)
    (nil))

(insn 56 55 57 2 ./CppStatUtilities.cc:426 (set (mem:DF (reg:DI 58 [ ivtmp.664 ]) [3 S8 A64])
        (reg:DF 89)) -1 (nil)
    (nil))

(note 57 56 58 2 ("./CppStatUtilities.cc") 425)

(insn 58 57 59 2 ./CppStatUtilities.cc:425 (parallel [
            (set (reg/v:SI 65 [ i ])
                (plus:SI (reg/v:SI 65 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 59 58 60 2 ./CppStatUtilities.cc:425 (parallel [
            (set (reg:DI 58 [ ivtmp.664 ])
                (plus:DI (reg:DI 58 [ ivtmp.664 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 60 59 61 2 ./CppStatUtilities.cc:425 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 65 [ i ])
            (reg/v:SI 77 [ nBins ]))) -1 (nil)
    (nil))

(jump_insn 61 60 62 2 ./CppStatUtilities.cc:425 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 50)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 62 61 63 3 451 "" [1 uses])

(note 63 62 64 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 64 63 65 3 ("./CppStatUtilities.cc") 429)

(insn 65 64 66 3 ./CppStatUtilities.cc:429 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 76 [ xSize ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 66 65 68 3 ./CppStatUtilities.cc:429 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 93)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 68 66 69 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 69 68 70 4 ./CppStatUtilities.cc:429 (set (reg:DI 59 [ ivtmp.659 ])
        (reg/v/f:DI 75 [ x ])) -1 (nil)
    (nil))

(insn 70 69 71 4 ./CppStatUtilities.cc:429 (set (reg/v:SI 64 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 71 70 72 5 456 "" [1 uses])

(note 72 71 73 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 73 72 74 5 ("./CppStatUtilities.cc") 430)

(insn 74 73 75 5 ./CppStatUtilities.cc:430 (set (reg:DF 91)
        (mem:DF (reg:DI 59 [ ivtmp.659 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 75 74 76 5 ./CppStatUtilities.cc:430 (set (reg:DF 90)
        (minus:DF (reg:DF 91)
            (reg/v:DF 71 [ xmin ]))) -1 (nil)
    (nil))

(insn 76 75 77 5 ./CppStatUtilities.cc:430 (set (reg:DF 92)
        (div:DF (reg:DF 90)
            (reg/v:DF 67 [ dx ]))) -1 (nil)
    (nil))

(insn 77 76 78 5 ./CppStatUtilities.cc:430 (set (reg:DF 21 xmm0)
        (reg:DF 92)) -1 (nil)
    (nil))

(call_insn/u 78 77 79 5 ./CppStatUtilities.cc:430 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("floor") [flags 0x41] <function_decl 0x2b5eb3c39300 floor>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 79 78 80 5 ./CppStatUtilities.cc:430 (set (reg:DF 73 [ D.35057 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 80 79 81 5 ./CppStatUtilities.cc:430 (set (reg:SI 93)
        (fix:SI (reg:DF 73 [ D.35057 ]))) -1 (nil)
    (nil))

(insn 81 80 82 5 ./CppStatUtilities.cc:430 (set (reg:DI 94)
        (sign_extend:DI (reg:SI 93))) -1 (nil)
    (nil))

(insn 82 81 83 5 ./CppStatUtilities.cc:430 (parallel [
            (set (reg:DI 95)
                (ashift:DI (reg:DI 94)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 83 82 84 5 ./CppStatUtilities.cc:430 (parallel [
            (set (reg:DI 72 [ D.35062 ])
                (plus:DI (reg:DI 95)
                    (reg/v/f:DI 68 [ hy ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 84 83 85 5 ./CppStatUtilities.cc:430 (set (reg:DF 96)
        (mem:DF (reg:DI 72 [ D.35062 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 85 84 86 5 ./CppStatUtilities.cc:430 (set (reg:DF 97)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC8") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 86 85 87 5 ./CppStatUtilities.cc:430 (set (reg:DF 98)
        (plus:DF (reg:DF 96)
            (reg:DF 97))) -1 (nil)
    (nil))

(insn 87 86 88 5 ./CppStatUtilities.cc:430 (set (mem:DF (reg:DI 72 [ D.35062 ]) [3 S8 A64])
        (reg:DF 98)) -1 (nil)
    (nil))

(note 88 87 89 5 ("./CppStatUtilities.cc") 429)

(insn 89 88 90 5 ./CppStatUtilities.cc:429 (parallel [
            (set (reg/v:SI 64 [ i ])
                (plus:SI (reg/v:SI 64 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 90 89 91 5 ./CppStatUtilities.cc:429 (parallel [
            (set (reg:DI 59 [ ivtmp.659 ])
                (plus:DI (reg:DI 59 [ ivtmp.659 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 91 90 92 5 ./CppStatUtilities.cc:429 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 64 [ i ])
            (reg/v:SI 76 [ xSize ]))) -1 (nil)
    (nil))

(jump_insn 92 91 93 5 ./CppStatUtilities.cc:429 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(code_label 93 92 94 6 454 "" [1 uses])

(note 94 93 95 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 95 94 96 6 ("./CppStatUtilities.cc") 433)

(insn 96 95 97 6 ./CppStatUtilities.cc:433 (set (reg:DF 99)
        (float:DF (reg/v:SI 76 [ xSize ]))) -1 (nil)
    (nil))

(insn 97 96 98 6 ./CppStatUtilities.cc:433 (set (reg/v:DF 66 [ f ])
        (mult:DF (reg:DF 99)
            (reg/v:DF 67 [ dx ]))) -1 (nil)
    (nil))

(note 98 97 99 6 ("./CppStatUtilities.cc") 434)

(insn 99 98 100 6 ./CppStatUtilities.cc:434 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 77 [ nBins ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 100 99 102 6 ./CppStatUtilities.cc:434 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 117)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 102 100 103 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 103 102 104 7 ./CppStatUtilities.cc:434 (set (reg/v:SI 63 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 104 103 105 8 459 "" [1 uses])

(note 105 104 106 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 106 105 107 8 ("./CppStatUtilities.cc") 412)

(insn 107 106 108 8 ./CppStatUtilities.cc:412 (set (reg:DI 60 [ D.36601 ])
        (reg:DI 61 [ ivtmp.654 ])) -1 (nil)
    (nil))

(note 108 107 109 8 ("./CppStatUtilities.cc") 435)

(insn 109 108 110 8 ./CppStatUtilities.cc:435 (set (reg:DF 100)
        (mem:DF (reg:DI 60 [ D.36601 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 110 109 111 8 ./CppStatUtilities.cc:435 (set (reg:DF 101)
        (div:DF (reg:DF 100)
            (reg/v:DF 66 [ f ]))) -1 (nil)
    (nil))

(insn 111 110 112 8 ./CppStatUtilities.cc:435 (set (mem:DF (reg:DI 60 [ D.36601 ]) [3 S8 A64])
        (reg:DF 101)) -1 (nil)
    (nil))

(note 112 111 113 8 ("./CppStatUtilities.cc") 434)

(insn 113 112 114 8 ./CppStatUtilities.cc:434 (parallel [
            (set (reg/v:SI 63 [ i ])
                (plus:SI (reg/v:SI 63 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 114 113 115 8 ./CppStatUtilities.cc:434 (parallel [
            (set (reg:DI 61 [ ivtmp.654 ])
                (plus:DI (reg:DI 61 [ ivtmp.654 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 115 114 116 8 ./CppStatUtilities.cc:434 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 63 [ i ])
            (reg/v:SI 77 [ nBins ]))) -1 (nil)
    (nil))

(jump_insn 116 115 117 8 ./CppStatUtilities.cc:434 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 104)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(code_label 117 116 118 9 457 "" [1 uses])

(note 118 117 119 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 119 118 120 9 ("./CppStatUtilities.cc") 437)

(insn 120 119 121 9 ./CppStatUtilities.cc:437 (set (mem/f:DI (reg/v/f:DI 78 [ _hx ]) [22 S8 A64])
        (reg/v/f:DI 69 [ hx ])) -1 (nil)
    (nil))

(note 121 120 122 9 ("./CppStatUtilities.cc") 438)

(insn 122 121 123 9 ./CppStatUtilities.cc:438 (set (mem/f:DI (reg/v/f:DI 79 [ _hy ]) [22 S8 A64])
        (reg/v/f:DI 68 [ hy ])) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

(note 123 122 124 NOTE_INSN_FUNCTION_END)

(note 124 123 0 ("./CppStatUtilities.cc") 439)


;; Function void summaryStats_t::print() const (_ZNK14summaryStats_t5printEv)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("./CppStatUtilities.cc") 54)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 ./CppStatUtilities.cc:54 (set (reg/f:DI 80 [ this ])
        (reg:DI 5 di [ this ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("./CppStatUtilities.cc") 59)

(insn 9 8 10 0 ./CppStatUtilities.cc:59 (set (reg:DF 79 [ D.34446 ])
        (mem/s:DF (plus:DI (reg/f:DI 80 [ this ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./CppStatUtilities.cc:59 (set (reg:DF 78 [ D.34447 ])
        (mem/s:DF (plus:DI (reg/f:DI 80 [ this ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppStatUtilities.cc:59 (set (reg:DF 77 [ D.34448 ])
        (mem/s:DF (plus:DI (reg/f:DI 80 [ this ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppStatUtilities.cc:59 (set (reg:DF 76 [ D.34449 ])
        (mem/s:DF (plus:DI (reg/f:DI 80 [ this ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:59 (set (reg:DF 75 [ D.34450 ])
        (mem/s:DF (plus:DI (reg/f:DI 80 [ this ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppStatUtilities.cc:59 (set (reg:DF 74 [ D.34451 ])
        (mem/s:DF (plus:DI (reg/f:DI 80 [ this ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:59 (set (reg:DF 73 [ D.34452 ])
        (mem/s:DF (reg/f:DI 80 [ this ]) [3 <variable>.min+0 S8 A64])) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC9") [flags 0x2] <string_cst 0x2b5eb6a34840>)) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt4cerr") [flags 0x40] <var_decl 0x2b5eb61c5840 cerr>)) -1 (nil)
    (nil))

(call_insn 18 17 19 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 19 18 20 0 ./CppStatUtilities.cc:59 (set (reg:DI 72 [ D.34453 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 20 19 21 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 73 [ D.34452 ])) -1 (nil)
    (nil))

(insn 21 20 22 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 72 [ D.34453 ])) -1 (nil)
    (nil))

(call_insn 22 21 23 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 23 22 24 0 ./CppStatUtilities.cc:59 (set (reg:DI 71 [ D.34454 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 24 23 25 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC10") [flags 0x2] <string_cst 0x2b5eb6a34930>)) -1 (nil)
    (nil))

(insn 25 24 26 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 71 [ D.34454 ])) -1 (nil)
    (nil))

(call_insn 26 25 27 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 27 26 28 0 ./CppStatUtilities.cc:59 (set (reg:DI 70 [ D.34455 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 28 27 29 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 74 [ D.34451 ])) -1 (nil)
    (nil))

(insn 29 28 30 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 70 [ D.34455 ])) -1 (nil)
    (nil))

(call_insn 30 29 31 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 31 30 32 0 ./CppStatUtilities.cc:59 (set (reg:DI 69 [ D.34456 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 32 31 33 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC11") [flags 0x2] <string_cst 0x2b5eb6a34960>)) -1 (nil)
    (nil))

(insn 33 32 34 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 69 [ D.34456 ])) -1 (nil)
    (nil))

(call_insn 34 33 35 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 35 34 36 0 ./CppStatUtilities.cc:59 (set (reg:DI 68 [ D.34457 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 36 35 37 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 75 [ D.34450 ])) -1 (nil)
    (nil))

(insn 37 36 38 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 68 [ D.34457 ])) -1 (nil)
    (nil))

(call_insn 38 37 39 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 39 38 40 0 ./CppStatUtilities.cc:59 (set (reg:DI 67 [ D.34458 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 40 39 41 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC12") [flags 0x2] <string_cst 0x2b5eb6a34a20>)) -1 (nil)
    (nil))

(insn 41 40 42 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 67 [ D.34458 ])) -1 (nil)
    (nil))

(call_insn 42 41 43 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 43 42 44 0 ./CppStatUtilities.cc:59 (set (reg:DI 66 [ D.34459 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 44 43 45 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 76 [ D.34449 ])) -1 (nil)
    (nil))

(insn 45 44 46 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 66 [ D.34459 ])) -1 (nil)
    (nil))

(call_insn 46 45 47 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 47 46 48 0 ./CppStatUtilities.cc:59 (set (reg:DI 65 [ D.34460 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 48 47 49 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC13") [flags 0x2] <string_cst 0x2b5eb6a34ae0>)) -1 (nil)
    (nil))

(insn 49 48 50 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 65 [ D.34460 ])) -1 (nil)
    (nil))

(call_insn 50 49 51 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 51 50 52 0 ./CppStatUtilities.cc:59 (set (reg:DI 64 [ D.34461 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 52 51 53 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 77 [ D.34448 ])) -1 (nil)
    (nil))

(insn 53 52 54 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 64 [ D.34461 ])) -1 (nil)
    (nil))

(call_insn 54 53 55 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 55 54 56 0 ./CppStatUtilities.cc:59 (set (reg:DI 63 [ D.34462 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 56 55 57 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC14") [flags 0x2] <string_cst 0x2b5eb6a34b10>)) -1 (nil)
    (nil))

(insn 57 56 58 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 63 [ D.34462 ])) -1 (nil)
    (nil))

(call_insn 58 57 59 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 59 58 60 0 ./CppStatUtilities.cc:59 (set (reg:DI 62 [ D.34463 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 60 59 61 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 78 [ D.34447 ])) -1 (nil)
    (nil))

(insn 61 60 62 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 62 [ D.34463 ])) -1 (nil)
    (nil))

(call_insn 62 61 63 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 63 62 64 0 ./CppStatUtilities.cc:59 (set (reg:DI 61 [ D.34464 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 64 63 65 0 ./CppStatUtilities.cc:59 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC15") [flags 0x2] <string_cst 0x2b5eb6a34b40>)) -1 (nil)
    (nil))

(insn 65 64 66 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 61 [ D.34464 ])) -1 (nil)
    (nil))

(call_insn 66 65 67 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 67 66 68 0 ./CppStatUtilities.cc:59 (set (reg:DI 60 [ D.34465 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 68 67 69 0 ./CppStatUtilities.cc:59 (set (reg:DF 21 xmm0)
        (reg:DF 79 [ D.34446 ])) -1 (nil)
    (nil))

(insn 69 68 70 0 ./CppStatUtilities.cc:59 (set (reg:DI 5 di)
        (reg:DI 60 [ D.34465 ])) -1 (nil)
    (nil))

(call_insn 70 69 71 0 ./CppStatUtilities.cc:59 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 71 70 72 0 ./CppStatUtilities.cc:59 (set (reg:DI 58 [ this ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 72 71 73 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc") 72)

(insn 73 72 74 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 5 di)
        (reg:DI 58 [ this ])) -1 (nil)
    (nil))

(call_insn/j 74 73 75 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_") [flags 0x41] <function_decl 0x2b5eb6063d00 endl>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 0, registers live:
 (nil)

(barrier 75 74 77)

(note 77 75 78 NOTE_INSN_FUNCTION_END)

(note 78 77 0 ("./CppStatUtilities.cc") 60)


;; Function double* percentiles(const double*, int, const double*, int, int) (_Z11percentilesPKdiS0_ii)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Deleted label in block 4.
Forwarding edge 4->5 to 12 failed.
Deleted label in block 8.
Forwarding edge 8->9 to 13 failed.
Deleted label in block 15.
Deleted label in block 18.
Deleted label in block 20.
Forwarding edge 20->21 to 29 failed.
Forwarding edge 22->23 to 29 failed.
Forwarding edge 25->26 to 29 failed.
Deleted label in block 30.
Forwarding edge 30->31 to 38 failed.
Deleted label in block 35.
Merged 38 and 39 without moving.
Merged 38 and 40 without moving.


try_optimize_cfg iteration 2

Forwarding edge 4->5 to 12 failed.
Forwarding edge 8->9 to 13 failed.
Forwarding edge 20->21 to 29 failed.
Forwarding edge 22->23 to 29 failed.
Forwarding edge 25->26 to 29 failed.
Forwarding edge 30->31 to 38 failed.


try_optimize_cfg iteration 1

Forwarding edge 3->4 to 11 failed.
Forwarding edge 7->8 to 12 failed.
Forwarding edge 19->20 to 28 failed.
Forwarding edge 21->22 to 28 failed.
Forwarding edge 24->25 to 28 failed.
Forwarding edge 29->30 to 37 failed.
(note 1 0 11 ("./CppStatUtilities.cc") 347)

;; Start of basic block 0, registers live: (nil)
(note 11 1 4 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 4 11 5 0 ./CppStatUtilities.cc:347 (set (reg/v/f:DI 101 [ data ])
        (reg:DI 5 di [ data ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:347 (set (reg/v:SI 102 [ dataLen ])
        (reg:SI 4 si [ dataLen ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:347 (set (reg/v/f:DI 103 [ prob ])
        (reg:DI 1 dx [ prob ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:347 (set (reg/v:SI 104 [ probLen ])
        (reg:SI 2 cx [ probLen ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppStatUtilities.cc:347 (set (reg/v:SI 105 [ perLen ])
        (reg:SI 37 r8 [ perLen ])) -1 (nil)
    (nil))

(note 9 8 13 0 NOTE_INSN_FUNCTION_BEG)

(note 13 9 14 0 ("./CppStatUtilities.cc") 349)

(insn 14 13 15 0 ./CppStatUtilities.cc:349 (parallel [
            (set (reg:SI 106)
                (plus:SI (reg/v:SI 105 [ perLen ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:349 (set (reg:CC 17 flags)
        (compare:CC (reg:SI 106)
            (const_int 99 [0x63]))) -1 (nil)
    (nil))

(jump_insn 16 15 18 0 ./CppStatUtilities.cc:349 (set (pc)
        (if_then_else (gtu (reg:CC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 37)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 100 [0x64])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 18 16 19 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 19 18 20 1 ("./CppStatUtilities.cc") 356)

(insn 20 19 21 1 ./CppStatUtilities.cc:356 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 104 [ probLen ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 21 20 23 1 ./CppStatUtilities.cc:356 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 152)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 23 21 24 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 24 23 25 2 ("./CppStatUtilities.cc") 357)

(insn 25 24 26 2 ./CppStatUtilities.cc:357 (set (reg:DF 59 [ temp.768 ])
        (mem:DF (reg/v/f:DI 103 [ prob ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 26 25 27 2 ./CppStatUtilities.cc:357 (set (reg:DF 107)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC17") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(insn 27 26 28 2 ./CppStatUtilities.cc:357 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 107)
            (reg:DF 59 [ temp.768 ]))) -1 (nil)
    (nil))

(jump_insn 28 27 30 2 ./CppStatUtilities.cc:357 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 85)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 100 [0x64])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 30 28 31 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 31 30 32 3 ./CppStatUtilities.cc:357 (set (reg:DF 108)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC18") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 32 31 33 3 ./CppStatUtilities.cc:357 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 59 [ temp.768 ])
            (reg:DF 108))) -1 (nil)
    (nil))

(jump_insn 33 32 36 3 ./CppStatUtilities.cc:357 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 85)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 36 33 34 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(jump_insn 34 36 35 4 ./CppStatUtilities.cc:357 (set (pc)
        (label_ref 142)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

(barrier 35 34 37)

;; Start of basic block 5, registers live: (nil)
(code_label 37 35 38 5 548 "" [1 uses])

(note 38 37 39 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 39 38 40 5 ("./CppStatUtilities.cc") 352)

(insn 40 39 41 5 ./CppStatUtilities.cc:352 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC19") [flags 0x2] <string_cst 0x2b5eb6ab7840>)) -1 (nil)
    (nil))

(insn 41 40 42 5 ./CppStatUtilities.cc:352 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt4cerr") [flags 0x40] <var_decl 0x2b5eb61c5840 cerr>)) -1 (nil)
    (nil))

(call_insn 42 41 43 5 ./CppStatUtilities.cc:352 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 43 42 44 5 ./CppStatUtilities.cc:352 (set (reg:DI 99 [ D.34898 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 44 43 45 5 ./CppStatUtilities.cc:352 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 45 44 46 5 ./CppStatUtilities.cc:352 (set (reg:DI 5 di)
        (reg:DI 99 [ D.34898 ])) -1 (nil)
    (nil))

(call_insn 46 45 47 5 ./CppStatUtilities.cc:352 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 47 46 48 5 ./CppStatUtilities.cc:352 (set (reg:DI 98 [ D.34899 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 48 47 49 5 ./CppStatUtilities.cc:352 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC20") [flags 0x2] <string_cst 0x2b5eb6abb880>)) -1 (nil)
    (nil))

(insn 49 48 50 5 ./CppStatUtilities.cc:352 (set (reg:DI 5 di)
        (reg:DI 98 [ D.34899 ])) -1 (nil)
    (nil))

(call_insn 50 49 51 5 ./CppStatUtilities.cc:352 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 51 50 52 5 ./CppStatUtilities.cc:352 (set (reg:DI 97 [ D.34900 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 52 51 53 5 ./CppStatUtilities.cc:352 (set (reg:SI 4 si)
        (const_int 351 [0x15f])) -1 (nil)
    (nil))

(insn 53 52 54 5 ./CppStatUtilities.cc:352 (set (reg:DI 5 di)
        (reg:DI 97 [ D.34900 ])) -1 (nil)
    (nil))

(call_insn 54 53 55 5 ./CppStatUtilities.cc:352 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEi") [flags 0x41] <function_decl 0x2b5eb6055b00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 55 54 56 5 ./CppStatUtilities.cc:352 (set (reg:DI 96 [ D.34901 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 56 55 57 5 ./CppStatUtilities.cc:352 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC21") [flags 0x2] <string_cst 0x2b5eb6ab8f50>)) -1 (nil)
    (nil))

(insn 57 56 58 5 ./CppStatUtilities.cc:352 (set (reg:DI 5 di)
        (reg:DI 96 [ D.34901 ])) -1 (nil)
    (nil))

(call_insn 58 57 59 5 ./CppStatUtilities.cc:352 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 59 58 60 5 ./CppStatUtilities.cc:352 (set (reg:DI 72 [ this ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 60 59 61 5 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc") 72)

(insn 61 60 62 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 5 di)
        (reg:DI 72 [ this ])) -1 (nil)
    (nil))

(call_insn 62 61 63 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_") [flags 0x41] <function_decl 0x2b5eb6063d00 endl>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 63 62 64 5 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 73 [ D.35471 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 64 63 65 5 ("./CppStatUtilities.cc") 353)

(insn 65 64 66 5 ./CppStatUtilities.cc:353 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 66 65 67 5 ./CppStatUtilities.cc:353 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(barrier 67 66 68)

;; Start of basic block 6, registers live: (nil)
(code_label 68 67 69 6 555 "" [1 uses])

(note 69 68 70 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(insn 70 69 71 6 ./CppStatUtilities.cc:353 (set (reg/v/f:DI 61 [ prob.766 ])
        (reg:DI 62 [ ivtmp.763 ])) -1 (nil)
    (nil))

(note 71 70 72 6 ("./CppStatUtilities.cc") 357)

(insn 72 71 73 6 ./CppStatUtilities.cc:357 (set (reg:DF 95 [ D.34913 ])
        (mem:DF (reg/v/f:DI 61 [ prob.766 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 73 72 74 6 ./CppStatUtilities.cc:357 (set (reg:DF 109)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC17") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(insn 74 73 75 6 ./CppStatUtilities.cc:357 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 109)
            (reg:DF 95 [ D.34913 ]))) -1 (nil)
    (nil))

(jump_insn 75 74 77 6 ./CppStatUtilities.cc:357 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 89)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 100 [0x64])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 77 75 78 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 78 77 79 7 ./CppStatUtilities.cc:357 (parallel [
            (set (reg:DI 62 [ ivtmp.763 ])
                (plus:DI (reg:DI 62 [ ivtmp.763 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 79 78 80 7 ./CppStatUtilities.cc:357 (set (reg:DF 110)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC18") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 80 79 81 7 ./CppStatUtilities.cc:357 (set (reg:CCFPU 17 flags)
        (compare:CCFPU (reg:DF 95 [ D.34913 ])
            (reg:DF 110))) -1 (nil)
    (nil))

(jump_insn 81 80 84 7 ./CppStatUtilities.cc:357 (set (pc)
        (if_then_else (gt (reg:CCFPU 17 flags)
                (const_int 0 [0x0]))
            (label_ref 89)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 333 [0x14d])
        (nil)))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(note 84 81 82 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(jump_insn 82 84 83 8 ./CppStatUtilities.cc:357 (set (pc)
        (label_ref 146)) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

(barrier 83 82 85)

;; Start of basic block 9, registers live: (nil)
(code_label 85 83 86 9 552 "" [2 uses])

(note 86 85 87 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(insn 87 86 88 9 ./CppStatUtilities.cc:357 (set (reg/v/f:DI 61 [ prob.766 ])
        (reg/v/f:DI 103 [ prob ])) -1 (nil)
    (nil))

(insn 88 87 89 9 ./CppStatUtilities.cc:357 (set (reg/v:SI 79 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 89 88 90 10 556 "" [2 uses])

(note 90 89 91 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 91 90 92 10 ("./CppStatUtilities.cc") 360)

(insn 92 91 93 10 ./CppStatUtilities.cc:360 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC19") [flags 0x2] <string_cst 0x2b5eb6ab7840>)) -1 (nil)
    (nil))

(insn 93 92 94 10 ./CppStatUtilities.cc:360 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt4cerr") [flags 0x40] <var_decl 0x2b5eb61c5840 cerr>)) -1 (nil)
    (nil))

(call_insn 94 93 95 10 ./CppStatUtilities.cc:360 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 95 94 96 10 ./CppStatUtilities.cc:360 (set (reg:DI 94 [ D.34914 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 96 95 97 10 ./CppStatUtilities.cc:360 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 97 96 98 10 ./CppStatUtilities.cc:360 (set (reg:DI 5 di)
        (reg:DI 94 [ D.34914 ])) -1 (nil)
    (nil))

(call_insn 98 97 99 10 ./CppStatUtilities.cc:360 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 99 98 100 10 ./CppStatUtilities.cc:360 (set (reg:DI 93 [ D.34915 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 100 99 101 10 ./CppStatUtilities.cc:360 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC20") [flags 0x2] <string_cst 0x2b5eb6abb880>)) -1 (nil)
    (nil))

(insn 101 100 102 10 ./CppStatUtilities.cc:360 (set (reg:DI 5 di)
        (reg:DI 93 [ D.34915 ])) -1 (nil)
    (nil))

(call_insn 102 101 103 10 ./CppStatUtilities.cc:360 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 103 102 104 10 ./CppStatUtilities.cc:360 (set (reg:DI 92 [ D.34916 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 104 103 105 10 ./CppStatUtilities.cc:360 (set (reg:SI 4 si)
        (const_int 359 [0x167])) -1 (nil)
    (nil))

(insn 105 104 106 10 ./CppStatUtilities.cc:360 (set (reg:DI 5 di)
        (reg:DI 92 [ D.34916 ])) -1 (nil)
    (nil))

(call_insn 106 105 107 10 ./CppStatUtilities.cc:360 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEi") [flags 0x41] <function_decl 0x2b5eb6055b00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 107 106 108 10 ./CppStatUtilities.cc:360 (set (reg:DI 91 [ D.34917 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 108 107 109 10 ./CppStatUtilities.cc:360 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC22") [flags 0x2] <string_cst 0x2b5eb6539a10>)) -1 (nil)
    (nil))

(insn 109 108 110 10 ./CppStatUtilities.cc:360 (set (reg:DI 5 di)
        (reg:DI 91 [ D.34917 ])) -1 (nil)
    (nil))

(call_insn 110 109 111 10 ./CppStatUtilities.cc:360 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 111 110 112 10 ./CppStatUtilities.cc:360 (set (reg:DI 70 [ this ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 112 111 113 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc") 72)

(insn 113 112 114 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 5 di)
        (reg:DI 70 [ this ])) -1 (nil)
    (nil))

(call_insn 114 113 115 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_") [flags 0x41] <function_decl 0x2b5eb6063d00 endl>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 115 114 116 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 71 [ D.35476 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 116 115 117 10 ("./CppStatUtilities.cc") 361)

(insn 117 116 118 10 ./CppStatUtilities.cc:361 (set (reg:DF 60 [ temp.767 ])
        (mem:DF (reg/v/f:DI 61 [ prob.766 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 118 117 119 10 ./CppStatUtilities.cc:361 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC23") [flags 0x2] <string_cst 0x2b5eb6ab7b10>)) -1 (nil)
    (nil))

(insn 119 118 120 10 ./CppStatUtilities.cc:361 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt4cerr") [flags 0x40] <var_decl 0x2b5eb61c5840 cerr>)) -1 (nil)
    (nil))

(call_insn 120 119 121 10 ./CppStatUtilities.cc:361 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 121 120 122 10 ./CppStatUtilities.cc:361 (set (reg:DI 90 [ D.34919 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 122 121 123 10 ./CppStatUtilities.cc:361 (set (reg:SI 4 si)
        (reg/v:SI 79 [ i ])) -1 (nil)
    (nil))

(insn 123 122 124 10 ./CppStatUtilities.cc:361 (set (reg:DI 5 di)
        (reg:DI 90 [ D.34919 ])) -1 (nil)
    (nil))

(call_insn 124 123 125 10 ./CppStatUtilities.cc:361 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEi") [flags 0x41] <function_decl 0x2b5eb6055b00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 125 124 126 10 ./CppStatUtilities.cc:361 (set (reg:DI 89 [ D.34920 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 126 125 127 10 ./CppStatUtilities.cc:361 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC24") [flags 0x2] <string_cst 0x2b5eb6ab7b40>)) -1 (nil)
    (nil))

(insn 127 126 128 10 ./CppStatUtilities.cc:361 (set (reg:DI 5 di)
        (reg:DI 89 [ D.34920 ])) -1 (nil)
    (nil))

(call_insn 128 127 129 10 ./CppStatUtilities.cc:361 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZStlsISt11char_traitsIcEERSt13basic_ostreamIcT_ES5_PKc") [flags 0x41] <function_decl 0x2b5eb606a900 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (nil))))

(insn 129 128 130 10 ./CppStatUtilities.cc:361 (set (reg:DI 88 [ D.34921 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 130 129 131 10 ./CppStatUtilities.cc:361 (set (reg:DF 21 xmm0)
        (reg:DF 60 [ temp.767 ])) -1 (nil)
    (nil))

(insn 131 130 132 10 ./CppStatUtilities.cc:361 (set (reg:DI 5 di)
        (reg:DI 88 [ D.34921 ])) -1 (nil)
    (nil))

(call_insn 132 131 133 10 ./CppStatUtilities.cc:361 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZNSolsEd") [flags 0x41] <function_decl 0x2b5eb6055f00 operator<<>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
            (nil))))

(insn 133 132 134 10 ./CppStatUtilities.cc:361 (set (reg:DI 68 [ this ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 134 133 135 10 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc") 72)

(insn 135 134 136 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 5 di)
        (reg:DI 68 [ this ])) -1 (nil)
    (nil))

(call_insn 136 135 137 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_ZSt4endlIcSt11char_traitsIcEERSt13basic_ostreamIT_T0_ES6_") [flags 0x41] <function_decl 0x2b5eb6063d00 endl>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 137 136 138 10 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/bits/ostream.tcc:72 (set (reg:DI 69 [ D.35481 ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 138 137 139 10 ("./CppStatUtilities.cc") 362)

(insn 139 138 140 10 ./CppStatUtilities.cc:362 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 140 139 141 10 ./CppStatUtilities.cc:362 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 10, registers live:
 (nil)

(barrier 141 140 142)

;; Start of basic block 11, registers live: (nil)
(code_label 142 141 143 11 554 "" [1 uses])

(note 143 142 144 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 144 143 145 11 ./CppStatUtilities.cc:362 (parallel [
            (set (reg:DI 62 [ ivtmp.763 ])
                (plus:DI (reg/v/f:DI 103 [ prob ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 145 144 146 11 ./CppStatUtilities.cc:362 (set (reg/v:SI 79 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(code_label 146 145 147 12 558 "" [1 uses])

(note 147 146 148 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(note 148 147 149 12 ("./CppStatUtilities.cc") 356)

(insn 149 148 150 12 ./CppStatUtilities.cc:356 (parallel [
            (set (reg/v:SI 79 [ i ])
                (plus:SI (reg/v:SI 79 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 150 149 151 12 ./CppStatUtilities.cc:356 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 79 [ i ])
            (reg/v:SI 104 [ probLen ]))) -1 (nil)
    (nil))

(jump_insn 151 150 152 12 ./CppStatUtilities.cc:356 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 68)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9667 [0x25c3])
        (nil)))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(code_label 152 151 153 13 550 "" [1 uses])

(note 153 152 154 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(note 154 153 155 13 ("./CppStatUtilities.cc") 367)

(insn 155 154 156 13 ./CppStatUtilities.cc:367 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 105 [ perLen ])
            (const_int 99 [0x63]))) -1 (nil)
    (nil))

(jump_insn 156 155 158 13 ./CppStatUtilities.cc:367 (set (pc)
        (if_then_else (le (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 162)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 13, registers live:
 (nil)

;; Start of basic block 14, registers live: (nil)
(note 158 156 159 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(insn 159 158 160 14 ./CppStatUtilities.cc:367 (set (reg/v:SI 82 [ sampleSize ])
        (reg/v:SI 102 [ dataLen ])) -1 (nil)
    (nil))

(jump_insn 160 159 161 14 ./CppStatUtilities.cc:367 (set (pc)
        (label_ref 171)) -1 (nil)
    (nil))
;; End of basic block 14, registers live:
 (nil)

(barrier 161 160 162)

;; Start of basic block 15, registers live: (nil)
(code_label 162 161 163 15 559 "" [1 uses])

(note 163 162 164 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 164 163 165 15 ("./CppStatUtilities.cc") 368)

(insn 165 164 166 15 ./CppStatUtilities.cc:368 (set (reg:DF 111)
        (float:DF (reg/v:SI 105 [ perLen ]))) -1 (nil)
    (nil))

(insn 166 165 167 15 ./CppStatUtilities.cc:368 (set (reg:DF 112)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC25") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+2 [0x0.c8p+7])
        (nil)))

(insn 167 166 168 15 ./CppStatUtilities.cc:368 (set (reg:DF 113)
        (div:DF (reg:DF 111)
            (reg:DF 112))) -1 (nil)
    (nil))

(insn 168 167 169 15 ./CppStatUtilities.cc:368 (set (reg:DF 114)
        (float:DF (reg/v:SI 102 [ dataLen ]))) -1 (nil)
    (nil))

(insn 169 168 170 15 ./CppStatUtilities.cc:368 (set (reg:DF 115)
        (mult:DF (reg:DF 113)
            (reg:DF 114))) -1 (nil)
    (nil))

(insn 170 169 171 15 ./CppStatUtilities.cc:368 (set (reg/v:SI 82 [ sampleSize ])
        (fix:SI (reg:DF 115))) -1 (nil)
    (nil))
;; End of basic block 15, registers live:
 (nil)

;; Start of basic block 16, registers live: (nil)
(code_label 171 170 172 16 561 "" [1 uses])

(note 172 171 173 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(note 173 172 174 16 ("./CppStatUtilities.cc") 372)

(insn 174 173 175 16 ./CppStatUtilities.cc:372 (set (reg:DI 87 [ D.34927 ])
        (sign_extend:DI (reg/v:SI 82 [ sampleSize ]))) -1 (nil)
    (nil))

(insn 175 174 176 16 ./CppStatUtilities.cc:372 (parallel [
            (set (reg:DI 86 [ D.34928 ])
                (ashift:DI (reg:DI 87 [ D.34927 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 176 175 177 16 ./CppStatUtilities.cc:372 (set (reg:DI 5 di)
        (reg:DI 86 [ D.34928 ])) -1 (nil)
    (nil))

(call_insn 177 176 178 16 ./CppStatUtilities.cc:372 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 178 177 179 16 ./CppStatUtilities.cc:372 (set (reg/f:DI 116)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 116)
        (nil)))

(insn 179 178 180 16 ./CppStatUtilities.cc:372 (set (reg:DI 63 [ ivtmp.758 ])
        (reg/f:DI 116)) -1 (nil)
    (nil))

(note 180 179 181 16 ("./CppStatUtilities.cc") 373)

(insn 181 180 182 16 ./CppStatUtilities.cc:373 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 63 [ ivtmp.758 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 182 181 184 16 ./CppStatUtilities.cc:373 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 194)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(note 184 182 185 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(insn 185 184 186 17 ./CppStatUtilities.cc:373 (set (reg:SI 2 cx)
        (const_int 373 [0x175])) -1 (nil)
    (nil))

(insn 186 185 187 17 ./CppStatUtilities.cc:373 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 187 186 188 17 ./CppStatUtilities.cc:373 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 188 187 189 17 ./CppStatUtilities.cc:373 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 189 188 190 17 ./CppStatUtilities.cc:373 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 190 189 191 17 ./CppStatUtilities.cc:373 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 191 190 192 17 ./CppStatUtilities.cc:373 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 192 191 193 17 ./CppStatUtilities.cc:373 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 17, registers live:
 (nil)

(barrier 193 192 194)

;; Start of basic block 18, registers live: (nil)
(code_label 194 193 195 18 562 "" [1 uses])

(note 195 194 196 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(note 196 195 197 18 ("./CppStatUtilities.cc") 372)

(insn 197 196 198 18 ./CppStatUtilities.cc:372 (set (reg/v/f:DI 81 [ b ])
        (reg:DI 63 [ ivtmp.758 ])) -1 (nil)
    (nil))

(note 198 197 199 18 ("./CppStatUtilities.cc") 375)

(insn 199 198 200 18 ./CppStatUtilities.cc:375 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 105 [ perLen ])
            (const_int 100 [0x64]))) -1 (nil)
    (nil))

(jump_insn 200 199 202 18 ./CppStatUtilities.cc:375 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 209)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 18, registers live:
 (nil)

;; Start of basic block 19, registers live: (nil)
(note 202 200 203 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(note 203 202 204 19 ("./CppStatUtilities.cc") 377)

(insn 204 203 205 19 ./CppStatUtilities.cc:377 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 82 [ sampleSize ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 205 204 208 19 ./CppStatUtilities.cc:377 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 217)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 19, registers live:
 (nil)

;; Start of basic block 20, registers live: (nil)
(note 208 205 206 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(jump_insn 206 208 207 20 ./CppStatUtilities.cc:377 (set (pc)
        (label_ref 251)) -1 (nil)
    (nil))
;; End of basic block 20, registers live:
 (nil)

(barrier 207 206 209)

;; Start of basic block 21, registers live: (nil)
(code_label 209 207 210 21 564 "" [1 uses])

(note 210 209 211 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(note 211 210 212 21 ("./CppStatUtilities.cc") 383)

(insn 212 211 213 21 ./CppStatUtilities.cc:383 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 82 [ sampleSize ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 213 212 216 21 ./CppStatUtilities.cc:383 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 234)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 21, registers live:
 (nil)

;; Start of basic block 22, registers live: (nil)
(note 216 213 214 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(jump_insn 214 216 215 22 ./CppStatUtilities.cc:383 (set (pc)
        (label_ref 251)) -1 (nil)
    (nil))
;; End of basic block 22, registers live:
 (nil)

(barrier 215 214 217)

;; Start of basic block 23, registers live: (nil)
(code_label 217 215 218 23 566 "" [1 uses])

(note 218 217 219 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(insn 219 218 220 23 ./CppStatUtilities.cc:383 (set (reg/v:SI 78 [ k ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 220 219 221 23 ./CppStatUtilities.cc:383 (set (reg:DI 64 [ ivtmp.750 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 23, registers live:
 (nil)

;; Start of basic block 24, registers live: (nil)
(code_label 221 220 222 24 569 "" [1 uses])

(note 222 221 223 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(note 223 222 224 24 ("./CppStatUtilities.cc") 378)

(insn 224 223 225 24 ./CppStatUtilities.cc:378 (set (reg:DF 117)
        (mem:DF (plus:DI (mult:DI (reg:DI 64 [ ivtmp.750 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 101 [ data ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 225 224 226 24 ./CppStatUtilities.cc:378 (set (mem:DF (plus:DI (mult:DI (reg:DI 64 [ ivtmp.750 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 81 [ b ])) [3 S8 A64])
        (reg:DF 117)) -1 (nil)
    (nil))

(note 226 225 227 24 ("./CppStatUtilities.cc") 377)

(insn 227 226 228 24 ./CppStatUtilities.cc:377 (parallel [
            (set (reg/v:SI 78 [ k ])
                (plus:SI (reg/v:SI 78 [ k ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 228 227 229 24 ./CppStatUtilities.cc:377 (parallel [
            (set (reg:DI 64 [ ivtmp.750 ])
                (plus:DI (reg:DI 64 [ ivtmp.750 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 229 228 230 24 ./CppStatUtilities.cc:377 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 78 [ k ])
            (reg/v:SI 82 [ sampleSize ]))) -1 (nil)
    (nil))

(jump_insn 230 229 233 24 ./CppStatUtilities.cc:377 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 221)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(note 233 230 231 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(jump_insn 231 233 232 25 ./CppStatUtilities.cc:377 (set (pc)
        (label_ref 251)) -1 (nil)
    (nil))
;; End of basic block 25, registers live:
 (nil)

(barrier 232 231 234)

;; Start of basic block 26, registers live: (nil)
(code_label 234 232 235 26 568 "" [1 uses])

(note 235 234 236 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(insn 236 235 237 26 ./CppStatUtilities.cc:377 (set (reg/v:SI 77 [ k ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 26, registers live:
 (nil)

;; Start of basic block 27, registers live: (nil)
(code_label 237 236 238 27 570 "" [1 uses])

(note 238 237 239 27 [bb 27] NOTE_INSN_BASIC_BLOCK)

(note 239 238 240 27 ("./CppStatUtilities.cc") 384)

(call_insn 240 239 241 27 ./CppStatUtilities.cc:384 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("rand") [flags 0x41] <function_decl 0x2b5eb50eef00 rand>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (nil))

(insn 241 240 242 27 ./CppStatUtilities.cc:384 (set (reg:SI 85 [ D.34952 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 242 241 243 27 ./CppStatUtilities.cc:384 (parallel [
            (set (reg:SI 120)
                (div:SI (reg:SI 85 [ D.34952 ])
                    (reg/v:SI 102 [ dataLen ])))
            (set (reg:SI 119)
                (mod:SI (reg:SI 85 [ D.34952 ])
                    (reg/v:SI 102 [ dataLen ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 243 242 244 27 ./CppStatUtilities.cc:384 (set (reg:DI 121)
        (sign_extend:DI (reg:SI 119))) -1 (nil)
    (nil))

(insn 244 243 245 27 ./CppStatUtilities.cc:384 (set (reg:DF 122)
        (mem:DF (plus:DI (mult:DI (reg:DI 121)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 101 [ data ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 245 244 246 27 ./CppStatUtilities.cc:384 (set (mem:DF (reg:DI 63 [ ivtmp.758 ]) [3 S8 A64])
        (reg:DF 122)) -1 (nil)
    (nil))

(note 246 245 247 27 ("./CppStatUtilities.cc") 383)

(insn 247 246 248 27 ./CppStatUtilities.cc:383 (parallel [
            (set (reg/v:SI 77 [ k ])
                (plus:SI (reg/v:SI 77 [ k ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 248 247 249 27 ./CppStatUtilities.cc:383 (parallel [
            (set (reg:DI 63 [ ivtmp.758 ])
                (plus:DI (reg:DI 63 [ ivtmp.758 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 249 248 250 27 ./CppStatUtilities.cc:383 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 77 [ k ])
            (reg/v:SI 82 [ sampleSize ]))) -1 (nil)
    (nil))

(jump_insn 250 249 251 27 ./CppStatUtilities.cc:383 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 237)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 27, registers live:
 (nil)

;; Start of basic block 28, registers live: (nil)
(code_label 251 250 252 28 567 "" [3 uses])

(note 252 251 253 28 [bb 28] NOTE_INSN_BASIC_BLOCK)

(note 253 252 254 28 ("./CppStatUtilities.cc") 387)

(insn 254 253 255 28 ./CppStatUtilities.cc:387 (set (reg:DI 2 cx)
        (symbol_ref:DI ("compareDouble") [flags 0x41] <function_decl 0x2b5eb6409400 compareDouble>)) -1 (nil)
    (nil))

(insn 255 254 256 28 ./CppStatUtilities.cc:387 (set (reg:DI 1 dx)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(insn 256 255 257 28 ./CppStatUtilities.cc:387 (set (reg:DI 4 si)
        (reg:DI 87 [ D.34927 ])) -1 (nil)
    (nil))

(insn 257 256 258 28 ./CppStatUtilities.cc:387 (set (reg:DI 5 di)
        (reg/v/f:DI 81 [ b ])) -1 (nil)
    (nil))

(call_insn 258 257 259 28 ./CppStatUtilities.cc:387 (call (mem:QI (symbol_ref:DI ("qsort") [flags 0x41] <function_decl 0x2b5eb5107600 qsort>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))

(note 259 258 260 28 ("./CppStatUtilities.cc") 389)

(insn 260 259 261 28 ./CppStatUtilities.cc:389 (set (reg:DI 123)
        (sign_extend:DI (reg/v:SI 104 [ probLen ]))) -1 (nil)
    (nil))

(insn 261 260 262 28 ./CppStatUtilities.cc:389 (parallel [
            (set (reg:DI 124)
                (ashift:DI (reg:DI 123)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 262 261 263 28 ./CppStatUtilities.cc:389 (set (reg:DI 5 di)
        (reg:DI 124)) -1 (nil)
    (nil))

(call_insn 263 262 264 28 ./CppStatUtilities.cc:389 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 264 263 265 28 ./CppStatUtilities.cc:389 (set (reg/f:DI 125)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 125)
        (nil)))

(insn 265 264 266 28 ./CppStatUtilities.cc:389 (set (reg:DI 84 [ D.34962 ])
        (reg/f:DI 125)) -1 (nil)
    (nil))

(note 266 265 267 28 ("./CppStatUtilities.cc") 390)

(insn 267 266 268 28 ./CppStatUtilities.cc:390 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 84 [ D.34962 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 268 267 270 28 ./CppStatUtilities.cc:390 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 279)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 100 [0x64])
        (nil)))
;; End of basic block 28, registers live:
 (nil)

;; Start of basic block 29, registers live: (nil)
(note 270 268 271 29 [bb 29] NOTE_INSN_BASIC_BLOCK)

(note 271 270 272 29 ("./CppStatUtilities.cc") 389)

(insn 272 271 273 29 ./CppStatUtilities.cc:389 (set (reg/v/f:DI 80 [ perc ])
        (reg:DI 84 [ D.34962 ])) -1 (nil)
    (nil))

(note 273 272 274 29 ("./CppStatUtilities.cc") 392)

(insn 274 273 275 29 ./CppStatUtilities.cc:392 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 104 [ probLen ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 275 274 278 29 ./CppStatUtilities.cc:392 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 291)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 29, registers live:
 (nil)

;; Start of basic block 30, registers live: (nil)
(note 278 275 276 30 [bb 30] NOTE_INSN_BASIC_BLOCK)

(jump_insn 276 278 277 30 ./CppStatUtilities.cc:392 (set (pc)
        (label_ref 340)) -1 (nil)
    (nil))
;; End of basic block 30, registers live:
 (nil)

(barrier 277 276 279)

;; Start of basic block 31, registers live: (nil)
(code_label 279 277 280 31 571 "" [1 uses])

(note 280 279 281 31 [bb 31] NOTE_INSN_BASIC_BLOCK)

(note 281 280 282 31 ("./CppStatUtilities.cc") 390)

(insn 282 281 283 31 ./CppStatUtilities.cc:390 (set (reg:SI 2 cx)
        (const_int 390 [0x186])) -1 (nil)
    (nil))

(insn 283 282 284 31 ./CppStatUtilities.cc:390 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 284 283 285 31 ./CppStatUtilities.cc:390 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 285 284 286 31 ./CppStatUtilities.cc:390 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 286 285 287 31 ./CppStatUtilities.cc:390 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 287 286 288 31 ./CppStatUtilities.cc:390 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 288 287 289 31 ./CppStatUtilities.cc:390 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 289 288 290 31 ./CppStatUtilities.cc:390 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 31, registers live:
 (nil)

(barrier 290 289 291)

;; Start of basic block 32, registers live: (nil)
(code_label 291 290 292 32 573 "" [1 uses])

(note 292 291 293 32 [bb 32] NOTE_INSN_BASIC_BLOCK)

(insn 293 292 294 32 ./CppStatUtilities.cc:390 (parallel [
            (set (reg:SI 67 [ pretmp.736 ])
                (plus:SI (reg/v:SI 82 [ sampleSize ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 294 293 295 32 ./CppStatUtilities.cc:390 (set (reg:DF 66 [ pretmp.737 ])
        (float:DF (reg:SI 67 [ pretmp.736 ]))) -1 (nil)
    (nil))

(insn 295 294 296 32 ./CppStatUtilities.cc:390 (set (reg/v:SI 76 [ k ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 296 295 297 32 ./CppStatUtilities.cc:390 (set (reg:DI 65 [ ivtmp.744 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 32, registers live:
 (nil)

;; Start of basic block 33, registers live: (nil)
(code_label 297 296 298 33 575 "" [1 uses])

(note 298 297 299 33 [bb 33] NOTE_INSN_BASIC_BLOCK)

(note 299 298 300 33 ("./CppStatUtilities.cc") 394)

(insn 300 299 301 33 ./CppStatUtilities.cc:394 (set (reg/v:DF 75 [ x ])
        (mult:DF (reg:DF 66 [ pretmp.737 ])
            (mem:DF (plus:DI (mult:DI (reg:DI 65 [ ivtmp.744 ])
                        (const_int 8 [0x8]))
                    (reg/v/f:DI 103 [ prob ])) [3 S8 A64]))) -1 (nil)
    (nil))

(note 301 300 302 33 ("./CppStatUtilities.cc") 395)

(insn 302 301 303 33 ./CppStatUtilities.cc:395 (set (reg:DF 21 xmm0)
        (reg/v:DF 75 [ x ])) -1 (nil)
    (nil))

(call_insn/u 303 302 304 33 ./CppStatUtilities.cc:395 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("floor") [flags 0x41] <function_decl 0x2b5eb3c39300 floor>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DF 21 xmm0))
        (nil)))

(insn 304 303 305 33 ./CppStatUtilities.cc:395 (set (reg:DF 83 [ D.34976 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 305 304 306 33 ("./CppStatUtilities.cc") 397)

(insn 306 305 307 33 ./CppStatUtilities.cc:397 (set (reg:SI 126)
        (fix:SI (reg:DF 83 [ D.34976 ]))) -1 (nil)
    (nil))

(insn 307 306 308 33 ./CppStatUtilities.cc:397 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 126)
            (reg:SI 67 [ pretmp.736 ]))) -1 (nil)
    (nil))

(jump_insn 308 307 310 33 ./CppStatUtilities.cc:397 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 316)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 33, registers live:
 (nil)

;; Start of basic block 34, registers live: (nil)
(note 310 308 311 34 [bb 34] NOTE_INSN_BASIC_BLOCK)

(note 311 310 312 34 ("./CppStatUtilities.cc") 398)

(insn 312 311 313 34 ./CppStatUtilities.cc:398 (set (reg:DF 127)
        (mem:DF (plus:DI (plus:DI (reg/v/f:DI 81 [ b ])
                    (reg:DI 86 [ D.34928 ]))
                (const_int -8 [0xfffffffffffffff8])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 313 312 314 34 ./CppStatUtilities.cc:398 (set (mem:DF (plus:DI (mult:DI (reg:DI 65 [ ivtmp.744 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 80 [ perc ])) [3 S8 A64])
        (reg:DF 127)) -1 (nil)
    (nil))

(jump_insn 314 313 315 34 ./CppStatUtilities.cc:398 (set (pc)
        (label_ref 333)) -1 (nil)
    (nil))
;; End of basic block 34, registers live:
 (nil)

(barrier 315 314 316)

;; Start of basic block 35, registers live: (nil)
(code_label 316 315 317 35 576 "" [1 uses])

(note 317 316 318 35 [bb 35] NOTE_INSN_BASIC_BLOCK)

(note 318 317 319 35 ("./CppStatUtilities.cc") 395)

(insn 319 318 320 35 ./CppStatUtilities.cc:395 (set (reg/v:DF 74 [ t ])
        (minus:DF (reg/v:DF 75 [ x ])
            (reg:DF 83 [ D.34976 ]))) -1 (nil)
    (nil))

(note 320 319 321 35 ("./CppStatUtilities.cc") 400)

(insn 321 320 322 35 ./CppStatUtilities.cc:400 (set (reg:DF 129)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC18") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 322 321 323 35 ./CppStatUtilities.cc:400 (set (reg:DF 128)
        (minus:DF (reg:DF 129)
            (reg/v:DF 74 [ t ]))) -1 (nil)
    (nil))

(insn 323 322 324 35 ./CppStatUtilities.cc:400 (set (reg:SI 130)
        (fix:SI (reg/v:DF 75 [ x ]))) -1 (nil)
    (nil))

(insn 324 323 325 35 ./CppStatUtilities.cc:400 (set (reg:DI 131)
        (sign_extend:DI (reg:SI 130))) -1 (nil)
    (nil))

(insn 325 324 326 35 ./CppStatUtilities.cc:400 (set (reg:DF 132)
        (mult:DF (reg:DF 128)
            (mem:DF (plus:DI (mult:DI (reg:DI 131)
                        (const_int 8 [0x8]))
                    (reg/v/f:DI 81 [ b ])) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 326 325 327 35 ./CppStatUtilities.cc:400 (set (reg:DF 134)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC18") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+0 [0x0.8p+1])
        (nil)))

(insn 327 326 328 35 ./CppStatUtilities.cc:400 (set (reg:DF 133)
        (plus:DF (reg/v:DF 75 [ x ])
            (reg:DF 134))) -1 (nil)
    (nil))

(insn 328 327 329 35 ./CppStatUtilities.cc:400 (set (reg:SI 135)
        (fix:SI (reg:DF 133))) -1 (nil)
    (nil))

(insn 329 328 330 35 ./CppStatUtilities.cc:400 (set (reg:DI 136)
        (sign_extend:DI (reg:SI 135))) -1 (nil)
    (nil))

(insn 330 329 331 35 ./CppStatUtilities.cc:400 (set (reg:DF 137)
        (mult:DF (reg/v:DF 74 [ t ])
            (mem:DF (plus:DI (mult:DI (reg:DI 136)
                        (const_int 8 [0x8]))
                    (reg/v/f:DI 81 [ b ])) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 331 330 332 35 ./CppStatUtilities.cc:400 (set (reg:DF 138)
        (plus:DF (reg:DF 132)
            (reg:DF 137))) -1 (nil)
    (nil))

(insn 332 331 333 35 ./CppStatUtilities.cc:400 (set (mem:DF (plus:DI (mult:DI (reg:DI 65 [ ivtmp.744 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 80 [ perc ])) [3 S8 A64])
        (reg:DF 138)) -1 (nil)
    (nil))
;; End of basic block 35, registers live:
 (nil)

;; Start of basic block 36, registers live: (nil)
(code_label 333 332 334 36 578 "" [1 uses])

(note 334 333 335 36 [bb 36] NOTE_INSN_BASIC_BLOCK)

(note 335 334 336 36 ("./CppStatUtilities.cc") 392)

(insn 336 335 337 36 ./CppStatUtilities.cc:392 (parallel [
            (set (reg/v:SI 76 [ k ])
                (plus:SI (reg/v:SI 76 [ k ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 337 336 338 36 ./CppStatUtilities.cc:392 (parallel [
            (set (reg:DI 65 [ ivtmp.744 ])
                (plus:DI (reg:DI 65 [ ivtmp.744 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 338 337 339 36 ./CppStatUtilities.cc:392 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 76 [ k ])
            (reg/v:SI 104 [ probLen ]))) -1 (nil)
    (nil))

(jump_insn 339 338 340 36 ./CppStatUtilities.cc:392 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 297)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 36, registers live:
 (nil)

;; Start of basic block 37, registers live: (nil)
(code_label 340 339 341 37 574 "" [1 uses])

(note 341 340 342 37 [bb 37] NOTE_INSN_BASIC_BLOCK)

(note 342 341 343 37 ("./CppStatUtilities.cc") 403)

(insn 343 342 344 37 ./CppStatUtilities.cc:403 (set (reg:DI 5 di)
        (reg/v/f:DI 81 [ b ])) -1 (nil)
    (nil))

(call_insn 344 343 345 37 ./CppStatUtilities.cc:403 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 345 344 348 37 ./CppStatUtilities.cc:403 (set (reg:DI 100 [ <result> ])
        (reg/v/f:DI 80 [ perc ])) -1 (nil)
    (nil))

(note 348 345 349 37 NOTE_INSN_FUNCTION_END)

(note 349 348 351 37 ("./CppStatUtilities.cc") 406)

(insn 351 349 357 37 ./CppStatUtilities.cc:406 (set (reg/i:DI 0 ax)
        (reg:DI 100 [ <result> ])) -1 (nil)
    (nil))

(insn 357 351 0 37 ./CppStatUtilities.cc:406 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 37, registers live:
 (nil)


;; Function double median(double*, double*, int, int) (_Z6medianPdS_ii)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Deleted label in block 6.
Deleted label in block 8.
Deleted label in block 11.
Merged 14 and 15 without moving.
Merged 14 and 16 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 9 ("./CppStatUtilities.cc") 301)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 ./CppStatUtilities.cc:301 (set (reg/v/f:DI 71 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:301 (set (reg/v/f:DI 72 [ w ])
        (reg:DI 4 si [ w ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:301 (set (reg/v:SI 73 [ start ])
        (reg:SI 1 dx [ start ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:301 (set (reg/v:SI 74 [ end ])
        (reg:SI 2 cx [ end ])) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("./CppStatUtilities.cc") 304)

(insn 12 11 13 0 ./CppStatUtilities.cc:304 (parallel [
            (set (reg:SI 69 [ D.34825 ])
                (minus:SI (reg/v:SI 74 [ end ])
                    (reg/v:SI 73 [ start ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:304 (parallel [
            (set (reg/v:SI 63 [ len ])
                (plus:SI (reg:SI 69 [ D.34825 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 14 13 15 0 ("./CppStatUtilities.cc") 306)

(insn 15 14 16 0 ./CppStatUtilities.cc:306 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 63 [ len ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 16 15 18 0 ./CppStatUtilities.cc:306 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 24)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 18 16 19 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 19 18 20 1 ("./CppStatUtilities.cc") 307)

(insn 20 19 21 1 ./CppStatUtilities.cc:307 (set (reg:DI 75)
        (sign_extend:DI (reg/v:SI 73 [ start ]))) -1 (nil)
    (nil))

(insn 21 20 22 1 ./CppStatUtilities.cc:307 (set (reg/v:DF 64 [ median ])
        (mem:DF (plus:DI (mult:DI (reg:DI 75)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 71 [ a ])) [3 S8 A64])) -1 (nil)
    (nil))

(jump_insn 22 21 23 1 ./CppStatUtilities.cc:307 (set (pc)
        (label_ref 142)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 23 22 24)

;; Start of basic block 2, registers live: (nil)
(code_label 24 23 25 2 713 "" [1 uses])

(note 25 24 26 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 26 25 27 2 ("./CppStatUtilities.cc") 308)

(insn 27 26 28 2 ./CppStatUtilities.cc:308 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 63 [ len ])
            (const_int 2 [0x2]))) -1 (nil)
    (nil))

(jump_insn 28 27 30 2 ./CppStatUtilities.cc:308 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 40)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5120 [0x1400])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 30 28 31 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 31 30 32 3 ("./CppStatUtilities.cc") 309)

(insn 32 31 33 3 ./CppStatUtilities.cc:309 (set (reg:DI 76)
        (sign_extend:DI (reg/v:SI 73 [ start ]))) -1 (nil)
    (nil))

(insn 33 32 34 3 ./CppStatUtilities.cc:309 (set (reg:DI 77)
        (sign_extend:DI (reg/v:SI 74 [ end ]))) -1 (nil)
    (nil))

(insn 34 33 35 3 ./CppStatUtilities.cc:309 (set (reg:DF 79)
        (mem:DF (plus:DI (mult:DI (reg:DI 76)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 71 [ a ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 35 34 36 3 ./CppStatUtilities.cc:309 (set (reg:DF 78)
        (plus:DF (reg:DF 79)
            (mem:DF (plus:DI (mult:DI (reg:DI 77)
                        (const_int 8 [0x8]))
                    (reg/v/f:DI 71 [ a ])) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 36 35 37 3 ./CppStatUtilities.cc:309 (set (reg:DF 80)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC26") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 37 36 38 3 ./CppStatUtilities.cc:309 (set (reg/v:DF 64 [ median ])
        (mult:DF (reg:DF 78)
            (reg:DF 80))) -1 (nil)
    (nil))

(jump_insn 38 37 39 3 ./CppStatUtilities.cc:309 (set (pc)
        (label_ref 142)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 39 38 40)

;; Start of basic block 4, registers live: (nil)
(code_label 40 39 41 4 716 "" [1 uses])

(note 41 40 42 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 42 41 43 4 ("./CppStatUtilities.cc") 311)

(insn 43 42 44 4 ./CppStatUtilities.cc:311 (set (reg:DI 67 [ D.34838 ])
        (sign_extend:DI (reg/v:SI 63 [ len ]))) -1 (nil)
    (nil))

(insn 44 43 45 4 ./CppStatUtilities.cc:311 (parallel [
            (set (reg:DI 66 [ D.34839 ])
                (ashift:DI (reg:DI 67 [ D.34838 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 45 44 46 4 ./CppStatUtilities.cc:311 (set (reg:DI 5 di)
        (reg:DI 66 [ D.34839 ])) -1 (nil)
    (nil))

(call_insn 46 45 47 4 ./CppStatUtilities.cc:311 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 47 46 48 4 ./CppStatUtilities.cc:311 (set (reg/f:DI 81)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 81)
        (nil)))

(insn 48 47 49 4 ./CppStatUtilities.cc:311 (set (reg:DI 65 [ D.34840 ])
        (reg/f:DI 81)) -1 (nil)
    (nil))

(note 49 48 50 4 ("./CppStatUtilities.cc") 312)

(insn 50 49 51 4 ./CppStatUtilities.cc:312 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 65 [ D.34840 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 51 50 53 4 ./CppStatUtilities.cc:312 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 63)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 53 51 54 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 54 53 55 5 ./CppStatUtilities.cc:312 (set (reg:SI 2 cx)
        (const_int 312 [0x138])) -1 (nil)
    (nil))

(insn 55 54 56 5 ./CppStatUtilities.cc:312 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 56 55 57 5 ./CppStatUtilities.cc:312 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 57 56 58 5 ./CppStatUtilities.cc:312 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 58 57 59 5 ./CppStatUtilities.cc:312 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 59 58 60 5 ./CppStatUtilities.cc:312 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 60 59 61 5 ./CppStatUtilities.cc:312 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 61 60 62 5 ./CppStatUtilities.cc:312 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(barrier 62 61 63)

;; Start of basic block 6, registers live: (nil)
(code_label 63 62 64 6 718 "" [1 uses])

(note 64 63 65 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 65 64 66 6 ("./CppStatUtilities.cc") 311)

(insn 66 65 67 6 ./CppStatUtilities.cc:311 (set (reg/v/f:DI 62 [ b ])
        (reg:DI 65 [ D.34840 ])) -1 (nil)
    (nil))

(note 67 66 68 6 ("./CppStatUtilities.cc") 314)

(insn 68 67 69 6 ./CppStatUtilities.cc:314 (set (reg:DI 82)
        (sign_extend:DI (reg/v:SI 73 [ start ]))) -1 (nil)
    (nil))

(insn 69 68 70 6 ./CppStatUtilities.cc:314 (parallel [
            (set (reg:DI 68 [ D.34828 ])
                (ashift:DI (reg:DI 82)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 70 69 73 6 ./CppStatUtilities.cc:314 (parallel [
            (set (reg:DI 83)
                (plus:DI (reg:DI 68 [ D.34828 ])
                    (reg/v/f:DI 71 [ a ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 73 70 74 6 ./CppStatUtilities.cc:314 (set (reg:DI 86)
        (reg/v/f:DI 62 [ b ])) -1 (nil)
    (nil))

(insn 74 73 75 6 ./CppStatUtilities.cc:314 (set (reg:DI 87)
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 75 74 76 6 ./CppStatUtilities.cc:314 (set (reg:DI 88)
        (reg:DI 66 [ D.34839 ])) -1 (nil)
    (nil))

(insn 76 75 77 6 ./CppStatUtilities.cc:314 (set (reg:DI 1 dx)
        (reg:DI 88)) -1 (nil)
    (nil))

(insn 77 76 78 6 ./CppStatUtilities.cc:314 (set (reg:DI 4 si)
        (reg:DI 87)) -1 (nil)
    (nil))

(insn 78 77 79 6 ./CppStatUtilities.cc:314 (set (reg:DI 5 di)
        (reg:DI 86)) -1 (nil)
    (nil))

(call_insn 79 78 80 6 ./CppStatUtilities.cc:314 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memcpy") [flags 0x41] <function_decl 0x2b5eb6bd5500 memcpy>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 80 79 81 6 ./CppStatUtilities.cc:314 (set (reg:DI 89)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 81 80 82 6 ("./CppStatUtilities.cc") 316)

(insn 82 81 83 6 ./CppStatUtilities.cc:316 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 63 [ len ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 83 82 85 6 ./CppStatUtilities.cc:316 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 103)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 85 83 86 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 86 85 87 7 ./CppStatUtilities.cc:316 (parallel [
            (set (reg:DI 60 [ ivtmp.808 ])
                (plus:DI (reg:DI 68 [ D.34828 ])
                    (reg/v/f:DI 72 [ w ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 87 86 88 7 ./CppStatUtilities.cc:316 (set (reg:DI 61 [ ivtmp.804 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 88 87 89 8 722 "" [1 uses])

(note 89 88 90 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 90 89 91 8 ("./CppStatUtilities.cc") 301)

(insn 91 90 92 8 ./CppStatUtilities.cc:301 (set (reg:DI 59 [ D.37083 ])
        (reg:DI 61 [ ivtmp.804 ])) -1 (nil)
    (nil))

(note 92 91 93 8 ("./CppStatUtilities.cc") 317)

(insn 93 92 94 8 ./CppStatUtilities.cc:317 (set (reg:DF 90)
        (mem:DF (plus:DI (mult:DI (reg:DI 59 [ D.37083 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 62 [ b ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 94 93 95 8 ./CppStatUtilities.cc:317 (set (reg:DF 91)
        (mult:DF (reg:DF 90)
            (mem:DF (reg:DI 60 [ ivtmp.808 ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 95 94 96 8 ./CppStatUtilities.cc:317 (set (mem:DF (plus:DI (mult:DI (reg:DI 59 [ D.37083 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 62 [ b ])) [3 S8 A64])
        (reg:DF 91)) -1 (nil)
    (nil))

(insn 96 95 97 8 ./CppStatUtilities.cc:317 (parallel [
            (set (reg:DI 61 [ ivtmp.804 ])
                (plus:DI (reg:DI 61 [ ivtmp.804 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 97 96 98 8 ./CppStatUtilities.cc:317 (parallel [
            (set (reg:DI 60 [ ivtmp.808 ])
                (plus:DI (reg:DI 60 [ ivtmp.808 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 98 97 99 8 ("./CppStatUtilities.cc") 316)

(insn 99 98 100 8 ./CppStatUtilities.cc:316 (set (reg:DI 92)
        (zero_extend:DI (reg:SI 69 [ D.34825 ]))) -1 (nil)
    (nil))

(insn 100 99 101 8 ./CppStatUtilities.cc:316 (parallel [
            (set (reg:DI 93)
                (plus:DI (reg:DI 92)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 101 100 102 8 ./CppStatUtilities.cc:316 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 61 [ ivtmp.804 ])
            (reg:DI 93))) -1 (nil)
    (nil))

(jump_insn 102 101 103 8 ./CppStatUtilities.cc:316 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 88)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(code_label 103 102 104 9 720 "" [1 uses])

(note 104 103 105 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 105 104 106 9 ("./CppStatUtilities.cc") 319)

(insn 106 105 107 9 ./CppStatUtilities.cc:319 (set (reg:DI 2 cx)
        (symbol_ref:DI ("compareDouble") [flags 0x41] <function_decl 0x2b5eb6409400 compareDouble>)) -1 (nil)
    (nil))

(insn 107 106 108 9 ./CppStatUtilities.cc:319 (set (reg:DI 1 dx)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(insn 108 107 109 9 ./CppStatUtilities.cc:319 (set (reg:DI 4 si)
        (reg:DI 67 [ D.34838 ])) -1 (nil)
    (nil))

(insn 109 108 110 9 ./CppStatUtilities.cc:319 (set (reg:DI 5 di)
        (reg/v/f:DI 62 [ b ])) -1 (nil)
    (nil))

(call_insn 110 109 111 9 ./CppStatUtilities.cc:319 (call (mem:QI (symbol_ref:DI ("qsort") [flags 0x41] <function_decl 0x2b5eb5107600 qsort>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))

(note 111 110 112 9 ("./CppStatUtilities.cc") 321)

(insn 112 111 113 9 ./CppStatUtilities.cc:321 (parallel [
            (set (reg:SI 94)
                (and:SI (reg/v:SI 63 [ len ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 113 112 114 9 ./CppStatUtilities.cc:321 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 94)
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 114 113 116 9 ./CppStatUtilities.cc:321 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 129)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(note 116 114 117 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 117 116 118 10 ("./CppStatUtilities.cc") 322)

(insn 118 117 119 10 ./CppStatUtilities.cc:322 (parallel [
            (set (reg:SI 96)
                (lshiftrt:SI (reg/v:SI 63 [ len ])
                    (const_int 31 [0x1f])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 119 118 120 10 ./CppStatUtilities.cc:322 (parallel [
            (set (reg:SI 97)
                (plus:SI (reg:SI 96)
                    (reg/v:SI 63 [ len ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 120 119 121 10 ./CppStatUtilities.cc:322 (parallel [
            (set (reg:SI 98)
                (ashiftrt:SI (reg:SI 97)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:SI (reg/v:SI 63 [ len ])
            (const_int 2 [0x2]))
        (nil)))

(insn 121 120 122 10 ./CppStatUtilities.cc:322 (set (reg:DI 99)
        (sign_extend:DI (reg:SI 98))) -1 (nil)
    (nil))

(insn 122 121 123 10 ./CppStatUtilities.cc:322 (parallel [
            (set (reg:DI 58 [ temp.818 ])
                (ashift:DI (reg:DI 99)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 123 122 124 10 ./CppStatUtilities.cc:322 (set (reg:DF 101)
        (mem:DF (plus:DI (plus:DI (reg/v/f:DI 62 [ b ])
                    (reg:DI 58 [ temp.818 ]))
                (const_int -8 [0xfffffffffffffff8])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 124 123 125 10 ./CppStatUtilities.cc:322 (set (reg:DF 100)
        (plus:DF (reg:DF 101)
            (mem:DF (plus:DI (reg/v/f:DI 62 [ b ])
                    (reg:DI 58 [ temp.818 ])) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 125 124 126 10 ./CppStatUtilities.cc:322 (set (reg:DF 102)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC26") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 126 125 127 10 ./CppStatUtilities.cc:322 (set (reg/v:DF 64 [ median ])
        (mult:DF (reg:DF 100)
            (reg:DF 102))) -1 (nil)
    (nil))

(jump_insn 127 126 128 10 ./CppStatUtilities.cc:322 (set (pc)
        (label_ref 137)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)

(barrier 128 127 129)

;; Start of basic block 11, registers live: (nil)
(code_label 129 128 130 11 723 "" [1 uses])

(note 130 129 131 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(note 131 130 132 11 ("./CppStatUtilities.cc") 324)

(insn 132 131 133 11 ./CppStatUtilities.cc:324 (parallel [
            (set (reg:SI 104)
                (lshiftrt:SI (reg/v:SI 63 [ len ])
                    (const_int 31 [0x1f])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 133 132 134 11 ./CppStatUtilities.cc:324 (parallel [
            (set (reg:SI 105)
                (plus:SI (reg:SI 104)
                    (reg/v:SI 63 [ len ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 134 133 135 11 ./CppStatUtilities.cc:324 (parallel [
            (set (reg:SI 106)
                (ashiftrt:SI (reg:SI 105)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:SI (reg/v:SI 63 [ len ])
            (const_int 2 [0x2]))
        (nil)))

(insn 135 134 136 11 ./CppStatUtilities.cc:324 (set (reg:DI 107)
        (sign_extend:DI (reg:SI 106))) -1 (nil)
    (nil))

(insn 136 135 137 11 ./CppStatUtilities.cc:324 (set (reg/v:DF 64 [ median ])
        (mem:DF (plus:DI (mult:DI (reg:DI 107)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 62 [ b ])) [3 S8 A64])) -1 (nil)
    (nil))
;; End of basic block 11, registers live:
 (nil)

;; Start of basic block 12, registers live: (nil)
(code_label 137 136 138 12 725 "" [1 uses])

(note 138 137 139 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(note 139 138 140 12 ("./CppStatUtilities.cc") 326)

(insn 140 139 141 12 ./CppStatUtilities.cc:326 (set (reg:DI 5 di)
        (reg/v/f:DI 62 [ b ])) -1 (nil)
    (nil))

(call_insn 141 140 142 12 ./CppStatUtilities.cc:326 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(code_label 142 141 143 13 715 "" [2 uses])

(note 143 142 144 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(insn 144 143 147 13 ./CppStatUtilities.cc:326 (set (reg:DF 70 [ <result> ])
        (reg/v:DF 64 [ median ])) -1 (nil)
    (nil))

(note 147 144 148 13 NOTE_INSN_FUNCTION_END)

(note 148 147 150 13 ("./CppStatUtilities.cc") 329)

(insn 150 148 156 13 ./CppStatUtilities.cc:329 (set (reg/i:DF 21 xmm0)
        (reg:DF 70 [ <result> ])) -1 (nil)
    (nil))

(insn 156 150 0 13 ./CppStatUtilities.cc:329 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 13, registers live:
 (nil)


;; Function double median(double*, int, int) (_Z6medianPdii)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Deleted label in block 6.
Deleted label in block 8.
Merged 11 and 12 without moving.
Merged 11 and 13 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 8 ("./CppStatUtilities.cc") 271)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./CppStatUtilities.cc:271 (set (reg/v/f:DI 66 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:271 (set (reg/v:SI 67 [ start ])
        (reg:SI 4 si [ start ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:271 (set (reg/v:SI 68 [ end ])
        (reg:SI 1 dx [ end ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./CppStatUtilities.cc") 274)

(insn 11 10 12 0 ./CppStatUtilities.cc:274 (parallel [
            (set (reg:SI 69)
                (minus:SI (reg/v:SI 68 [ end ])
                    (reg/v:SI 67 [ start ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppStatUtilities.cc:274 (parallel [
            (set (reg/v:SI 60 [ len ])
                (plus:SI (reg:SI 69)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 13 12 14 0 ("./CppStatUtilities.cc") 276)

(insn 14 13 15 0 ./CppStatUtilities.cc:276 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 60 [ len ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 15 14 17 0 ./CppStatUtilities.cc:276 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 23)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 7100 [0x1bbc])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 17 15 18 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 18 17 19 1 ("./CppStatUtilities.cc") 277)

(insn 19 18 20 1 ./CppStatUtilities.cc:277 (set (reg:DI 70)
        (sign_extend:DI (reg/v:SI 67 [ start ]))) -1 (nil)
    (nil))

(insn 20 19 21 1 ./CppStatUtilities.cc:277 (set (reg/v:DF 61 [ median ])
        (mem:DF (plus:DI (mult:DI (reg:DI 70)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 66 [ a ])) [3 S8 A64])) -1 (nil)
    (nil))

(jump_insn 21 20 22 1 ./CppStatUtilities.cc:277 (set (pc)
        (label_ref 117)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 22 21 23)

;; Start of basic block 2, registers live: (nil)
(code_label 23 22 24 2 764 "" [1 uses])

(note 24 23 25 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 25 24 26 2 ("./CppStatUtilities.cc") 278)

(insn 26 25 27 2 ./CppStatUtilities.cc:278 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 60 [ len ])
            (const_int 2 [0x2]))) -1 (nil)
    (nil))

(jump_insn 27 26 29 2 ./CppStatUtilities.cc:278 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 39)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5120 [0x1400])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 29 27 30 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 30 29 31 3 ("./CppStatUtilities.cc") 279)

(insn 31 30 32 3 ./CppStatUtilities.cc:279 (set (reg:DI 71)
        (sign_extend:DI (reg/v:SI 67 [ start ]))) -1 (nil)
    (nil))

(insn 32 31 33 3 ./CppStatUtilities.cc:279 (set (reg:DI 72)
        (sign_extend:DI (reg/v:SI 68 [ end ]))) -1 (nil)
    (nil))

(insn 33 32 34 3 ./CppStatUtilities.cc:279 (set (reg:DF 74)
        (mem:DF (plus:DI (mult:DI (reg:DI 71)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 66 [ a ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 34 33 35 3 ./CppStatUtilities.cc:279 (set (reg:DF 73)
        (plus:DF (reg:DF 74)
            (mem:DF (plus:DI (mult:DI (reg:DI 72)
                        (const_int 8 [0x8]))
                    (reg/v/f:DI 66 [ a ])) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 35 34 36 3 ./CppStatUtilities.cc:279 (set (reg:DF 75)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC27") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 36 35 37 3 ./CppStatUtilities.cc:279 (set (reg/v:DF 61 [ median ])
        (mult:DF (reg:DF 73)
            (reg:DF 75))) -1 (nil)
    (nil))

(jump_insn 37 36 38 3 ./CppStatUtilities.cc:279 (set (pc)
        (label_ref 117)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 38 37 39)

;; Start of basic block 4, registers live: (nil)
(code_label 39 38 40 4 767 "" [1 uses])

(note 40 39 41 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 41 40 42 4 ("./CppStatUtilities.cc") 281)

(insn 42 41 43 4 ./CppStatUtilities.cc:281 (set (reg:DI 64 [ D.34796 ])
        (sign_extend:DI (reg/v:SI 60 [ len ]))) -1 (nil)
    (nil))

(insn 43 42 44 4 ./CppStatUtilities.cc:281 (parallel [
            (set (reg:DI 63 [ D.34797 ])
                (ashift:DI (reg:DI 64 [ D.34796 ])
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 44 43 45 4 ./CppStatUtilities.cc:281 (set (reg:DI 5 di)
        (reg:DI 63 [ D.34797 ])) -1 (nil)
    (nil))

(call_insn 45 44 46 4 ./CppStatUtilities.cc:281 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 46 45 47 4 ./CppStatUtilities.cc:281 (set (reg/f:DI 76)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 76)
        (nil)))

(insn 47 46 48 4 ./CppStatUtilities.cc:281 (set (reg:DI 62 [ D.34798 ])
        (reg/f:DI 76)) -1 (nil)
    (nil))

(note 48 47 49 4 ("./CppStatUtilities.cc") 282)

(insn 49 48 50 4 ./CppStatUtilities.cc:282 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 62 [ D.34798 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 50 49 52 4 ./CppStatUtilities.cc:282 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 62)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 52 50 53 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 53 52 54 5 ./CppStatUtilities.cc:282 (set (reg:SI 2 cx)
        (const_int 282 [0x11a])) -1 (nil)
    (nil))

(insn 54 53 55 5 ./CppStatUtilities.cc:282 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 55 54 56 5 ./CppStatUtilities.cc:282 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 56 55 57 5 ./CppStatUtilities.cc:282 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 57 56 58 5 ./CppStatUtilities.cc:282 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 58 57 59 5 ./CppStatUtilities.cc:282 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 59 58 60 5 ./CppStatUtilities.cc:282 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 60 59 61 5 ./CppStatUtilities.cc:282 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 5, registers live:
 (nil)

(barrier 61 60 62)

;; Start of basic block 6, registers live: (nil)
(code_label 62 61 63 6 769 "" [1 uses])

(note 63 62 64 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 64 63 65 6 ("./CppStatUtilities.cc") 281)

(insn 65 64 66 6 ./CppStatUtilities.cc:281 (set (reg/v/f:DI 59 [ b ])
        (reg:DI 62 [ D.34798 ])) -1 (nil)
    (nil))

(note 66 65 67 6 ("./CppStatUtilities.cc") 284)

(insn 67 66 68 6 ./CppStatUtilities.cc:284 (set (reg:DI 77)
        (sign_extend:DI (reg/v:SI 67 [ start ]))) -1 (nil)
    (nil))

(insn 68 67 69 6 ./CppStatUtilities.cc:284 (parallel [
            (set (reg:DI 78)
                (ashift:DI (reg:DI 77)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 69 68 72 6 ./CppStatUtilities.cc:284 (parallel [
            (set (reg:DI 79)
                (plus:DI (reg:DI 78)
                    (reg/v/f:DI 66 [ a ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 72 69 73 6 ./CppStatUtilities.cc:284 (set (reg:DI 82)
        (reg/v/f:DI 59 [ b ])) -1 (nil)
    (nil))

(insn 73 72 74 6 ./CppStatUtilities.cc:284 (set (reg:DI 83)
        (reg:DI 79)) -1 (nil)
    (nil))

(insn 74 73 75 6 ./CppStatUtilities.cc:284 (set (reg:DI 84)
        (reg:DI 63 [ D.34797 ])) -1 (nil)
    (nil))

(insn 75 74 76 6 ./CppStatUtilities.cc:284 (set (reg:DI 1 dx)
        (reg:DI 84)) -1 (nil)
    (nil))

(insn 76 75 77 6 ./CppStatUtilities.cc:284 (set (reg:DI 4 si)
        (reg:DI 83)) -1 (nil)
    (nil))

(insn 77 76 78 6 ./CppStatUtilities.cc:284 (set (reg:DI 5 di)
        (reg:DI 82)) -1 (nil)
    (nil))

(call_insn 78 77 79 6 ./CppStatUtilities.cc:284 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("memcpy") [flags 0x41] <function_decl 0x2b5eb6bd5500 memcpy>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (nil)))))

(insn 79 78 80 6 ./CppStatUtilities.cc:284 (set (reg:DI 85)
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 80 79 81 6 ("./CppStatUtilities.cc") 286)

(insn 81 80 82 6 ./CppStatUtilities.cc:286 (set (reg:DI 2 cx)
        (symbol_ref:DI ("compareDouble") [flags 0x41] <function_decl 0x2b5eb6409400 compareDouble>)) -1 (nil)
    (nil))

(insn 82 81 83 6 ./CppStatUtilities.cc:286 (set (reg:DI 1 dx)
        (const_int 8 [0x8])) -1 (nil)
    (nil))

(insn 83 82 84 6 ./CppStatUtilities.cc:286 (set (reg:DI 4 si)
        (reg:DI 64 [ D.34796 ])) -1 (nil)
    (nil))

(insn 84 83 85 6 ./CppStatUtilities.cc:286 (set (reg:DI 5 di)
        (reg/v/f:DI 59 [ b ])) -1 (nil)
    (nil))

(call_insn 85 84 86 6 ./CppStatUtilities.cc:286 (call (mem:QI (symbol_ref:DI ("qsort") [flags 0x41] <function_decl 0x2b5eb5107600 qsort>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:DI 2 cx))
                    (nil))))))

(note 86 85 87 6 ("./CppStatUtilities.cc") 288)

(insn 87 86 88 6 ./CppStatUtilities.cc:288 (parallel [
            (set (reg:SI 86)
                (and:SI (reg/v:SI 60 [ len ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 88 87 89 6 ./CppStatUtilities.cc:288 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:SI 86)
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 89 88 91 6 ./CppStatUtilities.cc:288 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 104)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 6, registers live:
 (nil)

;; Start of basic block 7, registers live: (nil)
(note 91 89 92 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 92 91 93 7 ("./CppStatUtilities.cc") 289)

(insn 93 92 94 7 ./CppStatUtilities.cc:289 (parallel [
            (set (reg:SI 88)
                (lshiftrt:SI (reg/v:SI 60 [ len ])
                    (const_int 31 [0x1f])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 94 93 95 7 ./CppStatUtilities.cc:289 (parallel [
            (set (reg:SI 89)
                (plus:SI (reg:SI 88)
                    (reg/v:SI 60 [ len ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 95 94 96 7 ./CppStatUtilities.cc:289 (parallel [
            (set (reg:SI 90)
                (ashiftrt:SI (reg:SI 89)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:SI (reg/v:SI 60 [ len ])
            (const_int 2 [0x2]))
        (nil)))

(insn 96 95 97 7 ./CppStatUtilities.cc:289 (set (reg:DI 91)
        (sign_extend:DI (reg:SI 90))) -1 (nil)
    (nil))

(insn 97 96 98 7 ./CppStatUtilities.cc:289 (parallel [
            (set (reg:DI 58 [ temp.860 ])
                (ashift:DI (reg:DI 91)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 98 97 99 7 ./CppStatUtilities.cc:289 (set (reg:DF 93)
        (mem:DF (plus:DI (plus:DI (reg/v/f:DI 59 [ b ])
                    (reg:DI 58 [ temp.860 ]))
                (const_int -8 [0xfffffffffffffff8])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 99 98 100 7 ./CppStatUtilities.cc:289 (set (reg:DF 92)
        (plus:DF (reg:DF 93)
            (mem:DF (plus:DI (reg/v/f:DI 59 [ b ])
                    (reg:DI 58 [ temp.860 ])) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 100 99 101 7 ./CppStatUtilities.cc:289 (set (reg:DF 94)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC27") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 5.0e-1 [0x0.8p+0])
        (nil)))

(insn 101 100 102 7 ./CppStatUtilities.cc:289 (set (reg/v:DF 61 [ median ])
        (mult:DF (reg:DF 92)
            (reg:DF 94))) -1 (nil)
    (nil))

(jump_insn 102 101 103 7 ./CppStatUtilities.cc:289 (set (pc)
        (label_ref 112)) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

(barrier 103 102 104)

;; Start of basic block 8, registers live: (nil)
(code_label 104 103 105 8 771 "" [1 uses])

(note 105 104 106 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 106 105 107 8 ("./CppStatUtilities.cc") 291)

(insn 107 106 108 8 ./CppStatUtilities.cc:291 (parallel [
            (set (reg:SI 96)
                (lshiftrt:SI (reg/v:SI 60 [ len ])
                    (const_int 31 [0x1f])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 108 107 109 8 ./CppStatUtilities.cc:291 (parallel [
            (set (reg:SI 97)
                (plus:SI (reg:SI 96)
                    (reg/v:SI 60 [ len ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 109 108 110 8 ./CppStatUtilities.cc:291 (parallel [
            (set (reg:SI 98)
                (ashiftrt:SI (reg:SI 97)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (div:SI (reg/v:SI 60 [ len ])
            (const_int 2 [0x2]))
        (nil)))

(insn 110 109 111 8 ./CppStatUtilities.cc:291 (set (reg:DI 99)
        (sign_extend:DI (reg:SI 98))) -1 (nil)
    (nil))

(insn 111 110 112 8 ./CppStatUtilities.cc:291 (set (reg/v:DF 61 [ median ])
        (mem:DF (plus:DI (mult:DI (reg:DI 99)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 59 [ b ])) [3 S8 A64])) -1 (nil)
    (nil))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(code_label 112 111 113 9 773 "" [1 uses])

(note 113 112 114 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 114 113 115 9 ("./CppStatUtilities.cc") 293)

(insn 115 114 116 9 ./CppStatUtilities.cc:293 (set (reg:DI 5 di)
        (reg/v/f:DI 59 [ b ])) -1 (nil)
    (nil))

(call_insn 116 115 117 9 ./CppStatUtilities.cc:293 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 9, registers live:
 (nil)

;; Start of basic block 10, registers live: (nil)
(code_label 117 116 118 10 766 "" [2 uses])

(note 118 117 119 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(insn 119 118 122 10 ./CppStatUtilities.cc:293 (set (reg:DF 65 [ <result> ])
        (reg/v:DF 61 [ median ])) -1 (nil)
    (nil))

(note 122 119 123 10 NOTE_INSN_FUNCTION_END)

(note 123 122 125 10 ("./CppStatUtilities.cc") 296)

(insn 125 123 131 10 ./CppStatUtilities.cc:296 (set (reg/i:DF 21 xmm0)
        (reg:DF 65 [ <result> ])) -1 (nil)
    (nil))

(insn 131 125 0 10 ./CppStatUtilities.cc:296 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 10, registers live:
 (nil)


;; Function double median(double*, int) (_Z6medianPdi)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 7 ("./CppStatUtilities.cc") 334)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./CppStatUtilities.cc:334 (set (reg/v/f:DI 60 [ data ])
        (reg:DI 5 di [ data ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:334 (set (reg/v:SI 61 [ dataLen ])
        (reg:SI 4 si [ dataLen ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./CppStatUtilities.cc") 336)

(insn 10 9 11 0 ./CppStatUtilities.cc:336 (parallel [
            (set (reg:SI 62)
                (plus:SI (reg/v:SI 61 [ dataLen ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 11 10 12 0 ./CppStatUtilities.cc:336 (set (reg:SI 1 dx)
        (reg:SI 62)) -1 (nil)
    (nil))

(insn 12 11 13 0 ./CppStatUtilities.cc:336 (set (reg:SI 4 si)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:336 (set (reg:DI 5 di)
        (reg/v/f:DI 60 [ data ])) -1 (nil)
    (nil))

(call_insn/j 14 13 15 0 ./CppStatUtilities.cc:336 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("_Z6medianPdii") [flags 0x3] <function_decl 0x2b5eb6361e00 median>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:SI 1 dx))
                (nil)))))
;; End of basic block 0, registers live:
 (nil)

(barrier 15 14 17)

(note 17 15 18 NOTE_INSN_FUNCTION_END)

(note 18 17 0 ("./CppStatUtilities.cc") 337)


;; Function void madInPlace(double*, int, double*, double&, double&, double) (_Z10madInPlacePdiS_RdS0_d)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 4 and 5 without moving.
Merged 4 and 6 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 11 ("./CppStatUtilities.cc") 257)

;; Start of basic block 0, registers live: (nil)
(note 11 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 11 4 0 ./CppStatUtilities.cc:257 (set (reg/v/f:DI 63 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:257 (set (reg/v:SI 64 [ n ])
        (reg:SI 4 si [ n ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:257 (set (reg/v/f:DI 65 [ res ])
        (reg:DI 1 dx [ res ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:257 (set (reg/v/f:DI 66 [ median ])
        (reg:DI 2 cx [ median ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:257 (set (reg/v/f:DI 67 [ mad ])
        (reg:DI 37 r8 [ mad ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppStatUtilities.cc:257 (set (reg/v:DF 68 [ constant ])
        (reg:DF 21 xmm0 [ constant ])) -1 (nil)
    (nil))

(note 9 8 13 0 NOTE_INSN_FUNCTION_BEG)

(note 13 9 14 0 ("./CppStatUtilities.cc") 259)

(insn 14 13 15 0 ./CppStatUtilities.cc:259 (set (reg:SI 4 si)
        (reg/v:SI 64 [ n ])) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:259 (set (reg:DI 5 di)
        (reg/v/f:DI 63 [ a ])) -1 (nil)
    (nil))

(call_insn 16 15 17 0 ./CppStatUtilities.cc:259 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 17 16 18 0 ./CppStatUtilities.cc:259 (set (reg:DF 62 [ D.34759 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:259 (set (mem:DF (reg/v/f:DI 66 [ median ]) [3 S8 A64])
        (reg:DF 62 [ D.34759 ])) -1 (nil)
    (nil))

(note 19 18 20 0 ("./CppStatUtilities.cc") 262)

(insn 20 19 21 0 ./CppStatUtilities.cc:262 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 64 [ n ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 21 20 23 0 ./CppStatUtilities.cc:262 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 41)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 23 21 24 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 24 23 25 1 ./CppStatUtilities.cc:262 (set (reg/v:SI 60 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 25 24 26 1 ./CppStatUtilities.cc:262 (set (reg:DI 59 [ ivtmp.921 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 26 25 27 2 783 "" [1 uses])

(note 27 26 28 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 28 27 29 2 ("./CppStatUtilities.cc") 257)

(insn 29 28 30 2 ./CppStatUtilities.cc:257 (set (reg:DI 58 [ D.37304 ])
        (reg:DI 59 [ ivtmp.921 ])) -1 (nil)
    (nil))

(note 30 29 31 2 ("./CppStatUtilities.cc") 263)

(insn 31 30 32 2 ./CppStatUtilities.cc:263 (set (reg:DF 70)
        (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.37304 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 63 [ a ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 32 31 33 2 ./CppStatUtilities.cc:263 (set (reg:DF 69)
        (minus:DF (reg:DF 70)
            (mem:DF (reg/v/f:DI 66 [ median ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 33 32 34 2 ./CppStatUtilities.cc:263 (set (reg:V2DF 71)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC28") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 34 33 35 2 ./CppStatUtilities.cc:263 (parallel [
            (set (reg:DF 72)
                (abs:DF (reg:DF 69)))
            (use (reg:V2DF 71))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 35 34 36 2 ./CppStatUtilities.cc:263 (set (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.37304 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 65 [ res ])) [3 S8 A64])
        (reg:DF 72)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 69))
        (nil)))

(note 36 35 37 2 ("./CppStatUtilities.cc") 262)

(insn 37 36 38 2 ./CppStatUtilities.cc:262 (parallel [
            (set (reg/v:SI 60 [ i ])
                (plus:SI (reg/v:SI 60 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 38 37 39 2 ./CppStatUtilities.cc:262 (parallel [
            (set (reg:DI 59 [ ivtmp.921 ])
                (plus:DI (reg:DI 59 [ ivtmp.921 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 39 38 40 2 ./CppStatUtilities.cc:262 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 60 [ i ])
            (reg/v:SI 64 [ n ]))) -1 (nil)
    (nil))

(jump_insn 40 39 41 2 ./CppStatUtilities.cc:262 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 26)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 41 40 42 3 781 "" [1 uses])

(note 42 41 43 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 43 42 44 3 ("./CppStatUtilities.cc") 265)

(insn 44 43 45 3 ./CppStatUtilities.cc:265 (set (reg:SI 4 si)
        (reg/v:SI 64 [ n ])) -1 (nil)
    (nil))

(insn 45 44 46 3 ./CppStatUtilities.cc:265 (set (reg:DI 5 di)
        (reg/v/f:DI 65 [ res ])) -1 (nil)
    (nil))

(call_insn 46 45 47 3 ./CppStatUtilities.cc:265 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 47 46 48 3 ./CppStatUtilities.cc:265 (set (reg:DF 61 [ D.34773 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 48 47 49 3 ./CppStatUtilities.cc:265 (set (reg:DF 73)
        (mult:DF (reg:DF 61 [ D.34773 ])
            (reg/v:DF 68 [ constant ]))) -1 (nil)
    (nil))

(insn 49 48 50 3 ./CppStatUtilities.cc:265 (set (mem:DF (reg/v/f:DI 67 [ mad ]) [3 S8 A64])
        (reg:DF 73)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(note 50 49 51 NOTE_INSN_FUNCTION_END)

(note 51 50 0 ("./CppStatUtilities.cc") 266)


;; Function double madInPlace(double*, int, double*, double) (_Z10madInPlacePdiS_d)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 4 and 5 without moving.
Merged 4 and 6 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 9 ("./CppStatUtilities.cc") 239)

;; Start of basic block 0, registers live: (nil)
(note 9 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 9 4 0 ./CppStatUtilities.cc:239 (set (reg/v/f:DI 64 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:239 (set (reg/v:SI 65 [ n ])
        (reg:SI 4 si [ n ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:239 (set (reg/v/f:DI 66 [ res ])
        (reg:DI 1 dx [ res ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:239 (set (reg/v:DF 67 [ constant ])
        (reg:DF 21 xmm0 [ constant ])) -1 (nil)
    (nil))

(note 7 6 11 0 NOTE_INSN_FUNCTION_BEG)

(note 11 7 12 0 ("./CppStatUtilities.cc") 241)

(insn 12 11 13 0 ./CppStatUtilities.cc:241 (set (reg:SI 4 si)
        (reg/v:SI 65 [ n ])) -1 (nil)
    (nil))

(insn 13 12 14 0 ./CppStatUtilities.cc:241 (set (reg:DI 5 di)
        (reg/v/f:DI 64 [ a ])) -1 (nil)
    (nil))

(call_insn 14 13 15 0 ./CppStatUtilities.cc:241 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 15 14 16 0 ./CppStatUtilities.cc:241 (set (reg/v:DF 62 [ med ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 16 15 17 0 ("./CppStatUtilities.cc") 245)

(insn 17 16 18 0 ./CppStatUtilities.cc:245 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 65 [ n ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 18 17 20 0 ./CppStatUtilities.cc:245 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 38)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 20 18 21 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 21 20 22 1 ./CppStatUtilities.cc:245 (set (reg/v:SI 61 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 22 21 23 1 ./CppStatUtilities.cc:245 (set (reg:DI 59 [ ivtmp.957 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 23 22 24 2 821 "" [1 uses])

(note 24 23 25 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 25 24 26 2 ("./CppStatUtilities.cc") 239)

(insn 26 25 27 2 ./CppStatUtilities.cc:239 (set (reg:DI 58 [ D.37379 ])
        (reg:DI 59 [ ivtmp.957 ])) -1 (nil)
    (nil))

(note 27 26 28 2 ("./CppStatUtilities.cc") 246)

(insn 28 27 29 2 ./CppStatUtilities.cc:246 (set (reg:DF 69)
        (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.37379 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 64 [ a ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 29 28 30 2 ./CppStatUtilities.cc:246 (set (reg:DF 68)
        (minus:DF (reg:DF 69)
            (reg/v:DF 62 [ med ]))) -1 (nil)
    (nil))

(insn 30 29 31 2 ./CppStatUtilities.cc:246 (set (reg:V2DF 70)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC29") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 31 30 32 2 ./CppStatUtilities.cc:246 (parallel [
            (set (reg:DF 71)
                (abs:DF (reg:DF 68)))
            (use (reg:V2DF 70))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 32 31 33 2 ./CppStatUtilities.cc:246 (set (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.37379 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 66 [ res ])) [3 S8 A64])
        (reg:DF 71)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 68))
        (nil)))

(note 33 32 34 2 ("./CppStatUtilities.cc") 245)

(insn 34 33 35 2 ./CppStatUtilities.cc:245 (parallel [
            (set (reg/v:SI 61 [ i ])
                (plus:SI (reg/v:SI 61 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 35 34 36 2 ./CppStatUtilities.cc:245 (parallel [
            (set (reg:DI 59 [ ivtmp.957 ])
                (plus:DI (reg:DI 59 [ ivtmp.957 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 36 35 37 2 ./CppStatUtilities.cc:245 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 61 [ i ])
            (reg/v:SI 65 [ n ]))) -1 (nil)
    (nil))

(jump_insn 37 36 38 2 ./CppStatUtilities.cc:245 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 23)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 38 37 39 3 819 "" [1 uses])

(note 39 38 40 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 40 39 41 3 ("./CppStatUtilities.cc") 248)

(insn 41 40 42 3 ./CppStatUtilities.cc:248 (set (reg:SI 4 si)
        (reg/v:SI 65 [ n ])) -1 (nil)
    (nil))

(insn 42 41 43 3 ./CppStatUtilities.cc:248 (set (reg:DI 5 di)
        (reg/v/f:DI 66 [ res ])) -1 (nil)
    (nil))

(call_insn 43 42 44 3 ./CppStatUtilities.cc:248 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 44 43 45 3 ./CppStatUtilities.cc:248 (set (reg/v:DF 60 [ mad ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 45 44 46 3 ./CppStatUtilities.cc:248 (set (reg:DF 72)
        (mult:DF (reg/v:DF 60 [ mad ])
            (reg/v:DF 67 [ constant ]))) -1 (nil)
    (nil))

(insn 46 45 49 3 ./CppStatUtilities.cc:248 (set (reg:DF 63 [ <result> ])
        (reg:DF 72)) -1 (nil)
    (nil))

(note 49 46 50 3 NOTE_INSN_FUNCTION_END)

(note 50 49 52 3 ("./CppStatUtilities.cc") 251)

(insn 52 50 58 3 ./CppStatUtilities.cc:251 (set (reg/i:DF 21 xmm0)
        (reg:DF 63 [ <result> ])) -1 (nil)
    (nil))

(insn 58 52 0 3 ./CppStatUtilities.cc:251 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)


;; Function double mad(double*, double*, int, int, double&, double) (_Z3madPdS_iiRdd)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Edge 6->8 redirected to 9
Forwarding edge 6->7 to 10 failed.
Forwarding edge 6->7 to 10 failed.
Deleting block 8.
Merged 11 and 12 without moving.
Merged 11 and 13 without moving.


try_optimize_cfg iteration 2

Forwarding edge 6->7 to 10 failed.


try_optimize_cfg iteration 1

Forwarding edge 5->6 to 8 failed.
(note 1 0 12 ("./CppStatUtilities.cc") 215)

;; Start of basic block 0, registers live: (nil)
(note 12 1 4 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 4 12 5 0 ./CppStatUtilities.cc:215 (set (reg/v/f:DI 71 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:215 (set (reg/v/f:DI 72 [ w ])
        (reg:DI 4 si [ w ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:215 (set (reg/v:SI 73 [ start ])
        (reg:SI 1 dx [ start ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:215 (set (reg/v:SI 74 [ end ])
        (reg:SI 2 cx [ end ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppStatUtilities.cc:215 (set (reg/v/f:DI 75 [ medianVal ])
        (reg:DI 37 r8 [ medianVal ])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./CppStatUtilities.cc:215 (set (reg/v:DF 76 [ constant ])
        (reg:DF 21 xmm0 [ constant ])) -1 (nil)
    (nil))

(note 10 9 14 0 NOTE_INSN_FUNCTION_BEG)

(note 14 10 15 0 ("./CppStatUtilities.cc") 218)

(insn 15 14 16 0 ./CppStatUtilities.cc:218 (parallel [
            (set (reg:SI 77)
                (minus:SI (reg/v:SI 74 [ end ])
                    (reg/v:SI 73 [ start ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 16 15 17 0 ./CppStatUtilities.cc:218 (parallel [
            (set (reg:SI 78)
                (plus:SI (reg:SI 77)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:218 (set (reg:DI 79)
        (sign_extend:DI (reg:SI 78))) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:218 (parallel [
            (set (reg:DI 80)
                (ashift:DI (reg:DI 79)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 19 18 20 0 ./CppStatUtilities.cc:218 (set (reg:DI 5 di)
        (reg:DI 80)) -1 (nil)
    (nil))

(call_insn 20 19 21 0 ./CppStatUtilities.cc:218 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 0 ./CppStatUtilities.cc:218 (set (reg/f:DI 81)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 81)
        (nil)))

(insn 22 21 23 0 ./CppStatUtilities.cc:218 (set (reg:DI 69 [ D.34700 ])
        (reg/f:DI 81)) -1 (nil)
    (nil))

(note 23 22 24 0 ("./CppStatUtilities.cc") 219)

(insn 24 23 25 0 ./CppStatUtilities.cc:219 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 69 [ D.34700 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 25 24 27 0 ./CppStatUtilities.cc:219 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 37)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 27 25 28 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 28 27 29 1 ./CppStatUtilities.cc:219 (set (reg:SI 2 cx)
        (const_int 219 [0xdb])) -1 (nil)
    (nil))

(insn 29 28 30 1 ./CppStatUtilities.cc:219 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 30 29 31 1 ./CppStatUtilities.cc:219 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 31 30 32 1 ./CppStatUtilities.cc:219 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 32 31 33 1 ./CppStatUtilities.cc:219 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 33 32 34 1 ./CppStatUtilities.cc:219 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 34 33 35 1 ./CppStatUtilities.cc:219 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 35 34 36 1 ./CppStatUtilities.cc:219 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 1, registers live:
 (nil)

(barrier 36 35 37)

;; Start of basic block 2, registers live: (nil)
(code_label 37 36 38 2 857 "" [1 uses])

(note 38 37 39 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 39 38 40 2 ("./CppStatUtilities.cc") 218)

(insn 40 39 41 2 ./CppStatUtilities.cc:218 (set (reg/v/f:DI 66 [ res ])
        (reg:DI 69 [ D.34700 ])) -1 (nil)
    (nil))

(note 41 40 42 2 ("./CppStatUtilities.cc") 221)

(insn 42 41 43 2 ./CppStatUtilities.cc:221 (set (reg:SI 2 cx)
        (reg/v:SI 74 [ end ])) -1 (nil)
    (nil))

(insn 43 42 44 2 ./CppStatUtilities.cc:221 (set (reg:SI 1 dx)
        (reg/v:SI 73 [ start ])) -1 (nil)
    (nil))

(insn 44 43 45 2 ./CppStatUtilities.cc:221 (set (reg:DI 4 si)
        (reg/v/f:DI 72 [ w ])) -1 (nil)
    (nil))

(insn 45 44 46 2 ./CppStatUtilities.cc:221 (set (reg:DI 5 di)
        (reg/v/f:DI 71 [ a ])) -1 (nil)
    (nil))

(call_insn 46 45 47 2 ./CppStatUtilities.cc:221 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("_Z6medianPdS_ii") [flags 0x3] <function_decl 0x2b5eb6361d00 median>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:SI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                    (nil))))))

(insn 47 46 48 2 ./CppStatUtilities.cc:221 (set (reg:DF 68 [ D.34704 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 48 47 49 2 ./CppStatUtilities.cc:221 (set (mem:DF (reg/v/f:DI 75 [ medianVal ]) [3 S8 A64])
        (reg:DF 68 [ D.34704 ])) -1 (nil)
    (nil))

(note 49 48 50 2 ("./CppStatUtilities.cc") 224)

(insn 50 49 51 2 ./CppStatUtilities.cc:224 (parallel [
            (set (reg/v:SI 59 [ end.998 ])
                (plus:SI (reg/v:SI 74 [ end ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 51 50 52 2 ("./CppStatUtilities.cc") 225)

(insn 52 51 53 2 ./CppStatUtilities.cc:225 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 73 [ start ])
            (reg/v:SI 59 [ end.998 ]))) -1 (nil)
    (nil))

(jump_insn 53 52 55 2 ./CppStatUtilities.cc:225 (set (pc)
        (if_then_else (lt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 59)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 55 53 56 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 56 55 57 3 ./CppStatUtilities.cc:225 (set (reg/v:SI 65 [ k ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(jump_insn 57 56 58 3 ./CppStatUtilities.cc:225 (set (pc)
        (label_ref 95)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

(barrier 58 57 59)

;; Start of basic block 4, registers live: (nil)
(code_label 59 58 60 4 859 "" [1 uses])

(note 60 59 61 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 61 60 62 4 ("./CppStatUtilities.cc") 215)

(insn 62 61 63 4 ./CppStatUtilities.cc:215 (set (reg:DI 82)
        (sign_extend:DI (reg/v:SI 73 [ start ]))) -1 (nil)
    (nil))

(insn 63 62 64 4 ./CppStatUtilities.cc:215 (parallel [
            (set (reg:DI 61 [ D.37477 ])
                (ashift:DI (reg:DI 82)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 64 63 65 4 ./CppStatUtilities.cc:215 (parallel [
            (set (reg:DI 62 [ ivtmp.995 ])
                (plus:DI (reg:DI 61 [ D.37477 ])
                    (reg/v/f:DI 72 [ w ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 65 64 66 4 ./CppStatUtilities.cc:215 (parallel [
            (set (reg:DI 60 [ ivtmp.997 ])
                (plus:DI (reg/v/f:DI 71 [ a ])
                    (reg:DI 61 [ D.37477 ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 66 65 67 4 ./CppStatUtilities.cc:215 (set (reg/v:SI 63 [ i ])
        (reg/v:SI 73 [ start ])) -1 (nil)
    (nil))

(insn 67 66 68 4 ./CppStatUtilities.cc:215 (set (reg/v:SI 65 [ k ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 68 67 69 5 862 "" [1 uses])

(note 69 68 70 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 70 69 71 5 ("./CppStatUtilities.cc") 226)

(insn 71 70 72 5 ./CppStatUtilities.cc:226 (set (reg:DF 67 [ D.34713 ])
        (mem:DF (reg:DI 62 [ ivtmp.995 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 72 71 73 5 ./CppStatUtilities.cc:226 (set (reg:DF 83)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC30") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(jump_insn 73 72 120 5 ./CppStatUtilities.cc:226 (parallel [
            (set (pc)
                (if_then_else (ne (reg:DF 67 [ D.34713 ])
                        (reg:DF 83))
                    (label_ref:DI 77)
                    (pc)))
            (clobber (reg:CCFP 18 fpsr))
            (clobber (reg:CCFP 17 flags))
        ]) 536 {*fp_jcc_1_sse} (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 5, registers live:
 (nil)

;; Start of basic block 6, registers live: (nil)
(note 120 73 74 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(jump_insn 74 120 75 6 ./CppStatUtilities.cc:226 (set (pc)
        (label_ref 87)) -1 (nil)
    (nil))
;; End of basic block 6, registers live:
 (nil)

(barrier 75 74 77)

;; Start of basic block 7, registers live: (nil)
(code_label 77 75 78 7 865 "" [1 uses])

(note 78 77 79 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(note 79 78 80 7 ("./CppStatUtilities.cc") 227)

(insn 80 79 81 7 ./CppStatUtilities.cc:227 (set (reg:DI 84)
        (sign_extend:DI (reg/v:SI 65 [ k ]))) -1 (nil)
    (nil))

(insn 81 80 82 7 ./CppStatUtilities.cc:227 (set (reg:DF 85)
        (mult:DF (reg:DF 67 [ D.34713 ])
            (mem:DF (reg:DI 60 [ ivtmp.997 ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 82 81 83 7 ./CppStatUtilities.cc:227 (set (reg:DF 86)
        (minus:DF (reg:DF 85)
            (mem:DF (reg/v/f:DI 75 [ medianVal ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 83 82 84 7 ./CppStatUtilities.cc:227 (set (reg:V2DF 87)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC31") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 84 83 85 7 ./CppStatUtilities.cc:227 (parallel [
            (set (reg:DF 88)
                (abs:DF (reg:DF 86)))
            (use (reg:V2DF 87))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 85 84 86 7 ./CppStatUtilities.cc:227 (set (mem:DF (plus:DI (mult:DI (reg:DI 84)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 66 [ res ])) [3 S8 A64])
        (reg:DF 88)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 86))
        (nil)))

(insn 86 85 87 7 ./CppStatUtilities.cc:227 (parallel [
            (set (reg/v:SI 65 [ k ])
                (plus:SI (reg/v:SI 65 [ k ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 87 86 88 8 863 "" [1 uses])

(note 88 87 89 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 89 88 90 8 ("./CppStatUtilities.cc") 225)

(insn 90 89 91 8 ./CppStatUtilities.cc:225 (parallel [
            (set (reg/v:SI 63 [ i ])
                (plus:SI (reg/v:SI 63 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 91 90 92 8 ./CppStatUtilities.cc:225 (parallel [
            (set (reg:DI 62 [ ivtmp.995 ])
                (plus:DI (reg:DI 62 [ ivtmp.995 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 92 91 93 8 ./CppStatUtilities.cc:225 (parallel [
            (set (reg:DI 60 [ ivtmp.997 ])
                (plus:DI (reg:DI 60 [ ivtmp.997 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 93 92 94 8 ./CppStatUtilities.cc:225 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 59 [ end.998 ])
            (reg/v:SI 63 [ i ]))) -1 (nil)
    (nil))

(jump_insn 94 93 95 8 ./CppStatUtilities.cc:225 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 68)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(code_label 95 94 96 9 861 "" [1 uses])

(note 96 95 97 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 97 96 98 9 ("./CppStatUtilities.cc") 229)

(insn 98 97 99 9 ./CppStatUtilities.cc:229 (set (reg:SI 4 si)
        (reg/v:SI 65 [ k ])) -1 (nil)
    (nil))

(insn 99 98 100 9 ./CppStatUtilities.cc:229 (set (reg:DI 5 di)
        (reg/v/f:DI 66 [ res ])) -1 (nil)
    (nil))

(call_insn 100 99 101 9 ./CppStatUtilities.cc:229 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 101 100 102 9 ./CppStatUtilities.cc:229 (set (reg/v:DF 64 [ mad ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 102 101 103 9 ("./CppStatUtilities.cc") 231)

(insn 103 102 104 9 ./CppStatUtilities.cc:231 (set (reg:DI 5 di)
        (reg/v/f:DI 66 [ res ])) -1 (nil)
    (nil))

(call_insn 104 103 105 9 ./CppStatUtilities.cc:231 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 105 104 106 9 ./CppStatUtilities.cc:231 (set (reg:DF 89)
        (mult:DF (reg/v:DF 64 [ mad ])
            (reg/v:DF 76 [ constant ]))) -1 (nil)
    (nil))

(insn 106 105 109 9 ./CppStatUtilities.cc:231 (set (reg:DF 70 [ <result> ])
        (reg:DF 89)) -1 (nil)
    (nil))

(note 109 106 110 9 NOTE_INSN_FUNCTION_END)

(note 110 109 112 9 ("./CppStatUtilities.cc") 234)

(insn 112 110 118 9 ./CppStatUtilities.cc:234 (set (reg/i:DF 21 xmm0)
        (reg:DF 70 [ <result> ])) -1 (nil)
    (nil))

(insn 118 112 0 9 ./CppStatUtilities.cc:234 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)


;; Function double mad(double*, int, int, double&, double) (_Z3madPdiiRdd)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Merged 6 and 7 without moving.
Merged 6 and 8 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 11 ("./CppStatUtilities.cc") 193)

;; Start of basic block 0, registers live: (nil)
(note 11 1 4 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 4 11 5 0 ./CppStatUtilities.cc:193 (set (reg/v/f:DI 68 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:193 (set (reg/v:SI 69 [ start ])
        (reg:SI 4 si [ start ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:193 (set (reg/v:SI 70 [ end ])
        (reg:SI 1 dx [ end ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:193 (set (reg/v/f:DI 71 [ medianVal ])
        (reg:DI 2 cx [ medianVal ])) -1 (nil)
    (nil))

(insn 8 7 9 0 ./CppStatUtilities.cc:193 (set (reg/v:DF 72 [ constant ])
        (reg:DF 21 xmm0 [ constant ])) -1 (nil)
    (nil))

(note 9 8 13 0 NOTE_INSN_FUNCTION_BEG)

(note 13 9 14 0 ("./CppStatUtilities.cc") 195)

(insn 14 13 15 0 ./CppStatUtilities.cc:195 (parallel [
            (set (reg:SI 73)
                (minus:SI (reg/v:SI 70 [ end ])
                    (reg/v:SI 69 [ start ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 15 14 16 0 ./CppStatUtilities.cc:195 (parallel [
            (set (reg/v:SI 64 [ len ])
                (plus:SI (reg:SI 73)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 16 15 17 0 ("./CppStatUtilities.cc") 196)

(insn 17 16 18 0 ./CppStatUtilities.cc:196 (set (reg:DI 74)
        (sign_extend:DI (reg/v:SI 64 [ len ]))) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:196 (parallel [
            (set (reg:DI 75)
                (ashift:DI (reg:DI 74)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 19 18 20 0 ./CppStatUtilities.cc:196 (set (reg:DI 5 di)
        (reg:DI 75)) -1 (nil)
    (nil))

(call_insn 20 19 21 0 ./CppStatUtilities.cc:196 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 0 ./CppStatUtilities.cc:196 (set (reg/f:DI 76)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 76)
        (nil)))

(insn 22 21 23 0 ./CppStatUtilities.cc:196 (set (reg:DI 66 [ D.34660 ])
        (reg/f:DI 76)) -1 (nil)
    (nil))

(note 23 22 24 0 ("./CppStatUtilities.cc") 197)

(insn 24 23 25 0 ./CppStatUtilities.cc:197 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 66 [ D.34660 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 25 24 27 0 ./CppStatUtilities.cc:197 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 37)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 27 25 28 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 28 27 29 1 ./CppStatUtilities.cc:197 (set (reg:SI 2 cx)
        (const_int 197 [0xc5])) -1 (nil)
    (nil))

(insn 29 28 30 1 ./CppStatUtilities.cc:197 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 30 29 31 1 ./CppStatUtilities.cc:197 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 31 30 32 1 ./CppStatUtilities.cc:197 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 32 31 33 1 ./CppStatUtilities.cc:197 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 33 32 34 1 ./CppStatUtilities.cc:197 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 34 33 35 1 ./CppStatUtilities.cc:197 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 35 34 36 1 ./CppStatUtilities.cc:197 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 1, registers live:
 (nil)

(barrier 36 35 37)

;; Start of basic block 2, registers live: (nil)
(code_label 37 36 38 2 948 "" [1 uses])

(note 38 37 39 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 39 38 40 2 ("./CppStatUtilities.cc") 196)

(insn 40 39 41 2 ./CppStatUtilities.cc:196 (set (reg/v/f:DI 63 [ res ])
        (reg:DI 66 [ D.34660 ])) -1 (nil)
    (nil))

(note 41 40 42 2 ("./CppStatUtilities.cc") 199)

(insn 42 41 43 2 ./CppStatUtilities.cc:199 (set (reg:SI 1 dx)
        (reg/v:SI 70 [ end ])) -1 (nil)
    (nil))

(insn 43 42 44 2 ./CppStatUtilities.cc:199 (set (reg:SI 4 si)
        (reg/v:SI 69 [ start ])) -1 (nil)
    (nil))

(insn 44 43 45 2 ./CppStatUtilities.cc:199 (set (reg:DI 5 di)
        (reg/v/f:DI 68 [ a ])) -1 (nil)
    (nil))

(call_insn 45 44 46 2 ./CppStatUtilities.cc:199 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("_Z6medianPdii") [flags 0x3] <function_decl 0x2b5eb6361e00 median>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:SI 1 dx))
                (nil)))))

(insn 46 45 47 2 ./CppStatUtilities.cc:199 (set (reg:DF 65 [ D.34664 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 47 46 48 2 ./CppStatUtilities.cc:199 (set (mem:DF (reg/v/f:DI 71 [ medianVal ]) [3 S8 A64])
        (reg:DF 65 [ D.34664 ])) -1 (nil)
    (nil))

(note 48 47 49 2 ("./CppStatUtilities.cc") 202)

(insn 49 48 50 2 ./CppStatUtilities.cc:202 (parallel [
            (set (reg/v:SI 59 [ end.1039 ])
                (plus:SI (reg/v:SI 70 [ end ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 50 49 51 2 ("./CppStatUtilities.cc") 203)

(insn 51 50 52 2 ./CppStatUtilities.cc:203 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 69 [ start ])
            (reg/v:SI 59 [ end.1039 ]))) -1 (nil)
    (nil))

(jump_insn 52 51 54 2 ./CppStatUtilities.cc:203 (set (pc)
        (if_then_else (ge (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 73)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 54 52 55 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 55 54 56 3 ./CppStatUtilities.cc:203 (set (reg:DI 77)
        (sign_extend:DI (reg/v:SI 69 [ start ]))) -1 (nil)
    (nil))

(insn 56 55 57 3 ./CppStatUtilities.cc:203 (parallel [
            (set (reg:DI 78)
                (ashift:DI (reg:DI 77)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 57 56 58 3 ./CppStatUtilities.cc:203 (parallel [
            (set (reg:DI 60 [ ivtmp.1036 ])
                (plus:DI (reg:DI 78)
                    (reg/v/f:DI 68 [ a ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 58 57 59 3 ./CppStatUtilities.cc:203 (set (reg:DI 61 [ ivtmp.1033 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 59 58 60 4 952 "" [1 uses])

(note 60 59 61 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 61 60 62 4 ("./CppStatUtilities.cc") 204)

(insn 62 61 63 4 ./CppStatUtilities.cc:204 (set (reg:DF 80)
        (mem:DF (reg:DI 60 [ ivtmp.1036 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 63 62 64 4 ./CppStatUtilities.cc:204 (set (reg:DF 79)
        (minus:DF (reg:DF 80)
            (mem:DF (reg/v/f:DI 71 [ medianVal ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 64 63 65 4 ./CppStatUtilities.cc:204 (set (reg:V2DF 81)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC32") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 65 64 66 4 ./CppStatUtilities.cc:204 (parallel [
            (set (reg:DF 82)
                (abs:DF (reg:DF 79)))
            (use (reg:V2DF 81))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 66 65 67 4 ./CppStatUtilities.cc:204 (set (mem:DF (plus:DI (mult:DI (reg:DI 61 [ ivtmp.1033 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 63 [ res ])) [3 S8 A64])
        (reg:DF 82)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 79))
        (nil)))

(insn 67 66 68 4 ./CppStatUtilities.cc:204 (parallel [
            (set (reg:DI 61 [ ivtmp.1033 ])
                (plus:DI (reg:DI 61 [ ivtmp.1033 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 68 67 69 4 ./CppStatUtilities.cc:204 (parallel [
            (set (reg:DI 60 [ ivtmp.1036 ])
                (plus:DI (reg:DI 60 [ ivtmp.1036 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 69 68 70 4 ("./CppStatUtilities.cc") 203)

(insn 70 69 71 4 ./CppStatUtilities.cc:203 (parallel [
            (set (reg:SI 83)
                (plus:SI (reg/v:SI 69 [ start ])
                    (subreg:SI (reg:DI 61 [ ivtmp.1033 ]) 0)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 71 70 72 4 ./CppStatUtilities.cc:203 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 59 [ end.1039 ])
            (reg:SI 83))) -1 (nil)
    (nil))

(jump_insn 72 71 73 4 ./CppStatUtilities.cc:203 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 59)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 73 72 74 5 950 "" [1 uses])

(note 74 73 75 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 75 74 76 5 ("./CppStatUtilities.cc") 206)

(insn 76 75 77 5 ./CppStatUtilities.cc:206 (set (reg:SI 4 si)
        (reg/v:SI 64 [ len ])) -1 (nil)
    (nil))

(insn 77 76 78 5 ./CppStatUtilities.cc:206 (set (reg:DI 5 di)
        (reg/v/f:DI 63 [ res ])) -1 (nil)
    (nil))

(call_insn 78 77 79 5 ./CppStatUtilities.cc:206 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 79 78 80 5 ./CppStatUtilities.cc:206 (set (reg/v:DF 62 [ mad ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 80 79 81 5 ("./CppStatUtilities.cc") 208)

(insn 81 80 82 5 ./CppStatUtilities.cc:208 (set (reg:DI 5 di)
        (reg/v/f:DI 63 [ res ])) -1 (nil)
    (nil))

(call_insn 82 81 83 5 ./CppStatUtilities.cc:208 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 83 82 84 5 ./CppStatUtilities.cc:208 (set (reg:DF 84)
        (mult:DF (reg/v:DF 62 [ mad ])
            (reg/v:DF 72 [ constant ]))) -1 (nil)
    (nil))

(insn 84 83 87 5 ./CppStatUtilities.cc:208 (set (reg:DF 67 [ <result> ])
        (reg:DF 84)) -1 (nil)
    (nil))

(note 87 84 88 5 NOTE_INSN_FUNCTION_END)

(note 88 87 90 5 ("./CppStatUtilities.cc") 211)

(insn 90 88 96 5 ./CppStatUtilities.cc:211 (set (reg/i:DF 21 xmm0)
        (reg:DF 67 [ <result> ])) -1 (nil)
    (nil))

(insn 96 90 0 5 ./CppStatUtilities.cc:211 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)


;; Function double mad(double*, int, int, double) (_Z3madPdiid)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 4.
Merged 6 and 7 without moving.
Merged 6 and 8 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 10 ("./CppStatUtilities.cc") 171)

;; Start of basic block 0, registers live: (nil)
(note 10 1 4 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 4 10 5 0 ./CppStatUtilities.cc:171 (set (reg/v/f:DI 68 [ a ])
        (reg:DI 5 di [ a ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:171 (set (reg/v:SI 69 [ start ])
        (reg:SI 4 si [ start ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:171 (set (reg/v:SI 70 [ end ])
        (reg:SI 1 dx [ end ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:171 (set (reg/v:DF 71 [ constant ])
        (reg:DF 21 xmm0 [ constant ])) -1 (nil)
    (nil))

(note 8 7 12 0 NOTE_INSN_FUNCTION_BEG)

(note 12 8 13 0 ("./CppStatUtilities.cc") 173)

(insn 13 12 14 0 ./CppStatUtilities.cc:173 (parallel [
            (set (reg:SI 72)
                (minus:SI (reg/v:SI 70 [ end ])
                    (reg/v:SI 69 [ start ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 14 13 15 0 ./CppStatUtilities.cc:173 (parallel [
            (set (reg/v:SI 65 [ len ])
                (plus:SI (reg:SI 72)
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 15 14 16 0 ("./CppStatUtilities.cc") 174)

(insn 16 15 17 0 ./CppStatUtilities.cc:174 (set (reg:DI 73)
        (sign_extend:DI (reg/v:SI 65 [ len ]))) -1 (nil)
    (nil))

(insn 17 16 18 0 ./CppStatUtilities.cc:174 (parallel [
            (set (reg:DI 74)
                (ashift:DI (reg:DI 73)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 18 17 19 0 ./CppStatUtilities.cc:174 (set (reg:DI 5 di)
        (reg:DI 74)) -1 (nil)
    (nil))

(call_insn 19 18 20 0 ./CppStatUtilities.cc:174 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 20 19 21 0 ./CppStatUtilities.cc:174 (set (reg/f:DI 75)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 75)
        (nil)))

(insn 21 20 22 0 ./CppStatUtilities.cc:174 (set (reg:DI 66 [ D.34623 ])
        (reg/f:DI 75)) -1 (nil)
    (nil))

(note 22 21 23 0 ("./CppStatUtilities.cc") 175)

(insn 23 22 24 0 ./CppStatUtilities.cc:175 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 66 [ D.34623 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 24 23 26 0 ./CppStatUtilities.cc:175 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 36)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 26 24 27 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 27 26 28 1 ./CppStatUtilities.cc:175 (set (reg:SI 2 cx)
        (const_int 175 [0xaf])) -1 (nil)
    (nil))

(insn 28 27 29 1 ./CppStatUtilities.cc:175 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 29 28 30 1 ./CppStatUtilities.cc:175 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 30 29 31 1 ./CppStatUtilities.cc:175 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 31 30 32 1 ./CppStatUtilities.cc:175 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 32 31 33 1 ./CppStatUtilities.cc:175 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 33 32 34 1 ./CppStatUtilities.cc:175 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 34 33 35 1 ./CppStatUtilities.cc:175 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 1, registers live:
 (nil)

(barrier 35 34 36)

;; Start of basic block 2, registers live: (nil)
(code_label 36 35 37 2 989 "" [1 uses])

(note 37 36 38 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 38 37 39 2 ("./CppStatUtilities.cc") 174)

(insn 39 38 40 2 ./CppStatUtilities.cc:174 (set (reg/v/f:DI 64 [ res ])
        (reg:DI 66 [ D.34623 ])) -1 (nil)
    (nil))

(note 40 39 41 2 ("./CppStatUtilities.cc") 177)

(insn 41 40 42 2 ./CppStatUtilities.cc:177 (set (reg:SI 1 dx)
        (reg/v:SI 70 [ end ])) -1 (nil)
    (nil))

(insn 42 41 43 2 ./CppStatUtilities.cc:177 (set (reg:SI 4 si)
        (reg/v:SI 69 [ start ])) -1 (nil)
    (nil))

(insn 43 42 44 2 ./CppStatUtilities.cc:177 (set (reg:DI 5 di)
        (reg/v/f:DI 68 [ a ])) -1 (nil)
    (nil))

(call_insn 44 43 45 2 ./CppStatUtilities.cc:177 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("_Z6medianPdii") [flags 0x3] <function_decl 0x2b5eb6361e00 median>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:SI 1 dx))
                (nil)))))

(insn 45 44 46 2 ./CppStatUtilities.cc:177 (set (reg/v:DF 63 [ med ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 46 45 47 2 ("./CppStatUtilities.cc") 180)

(insn 47 46 48 2 ./CppStatUtilities.cc:180 (parallel [
            (set (reg/v:SI 59 [ end.1080 ])
                (plus:SI (reg/v:SI 70 [ end ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 48 47 49 2 ("./CppStatUtilities.cc") 181)

(insn 49 48 50 2 ./CppStatUtilities.cc:181 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 69 [ start ])
            (reg/v:SI 59 [ end.1080 ]))) -1 (nil)
    (nil))

(jump_insn 50 49 52 2 ./CppStatUtilities.cc:181 (set (pc)
        (if_then_else (ge (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 71)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(note 52 50 53 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 53 52 54 3 ./CppStatUtilities.cc:181 (set (reg:DI 76)
        (sign_extend:DI (reg/v:SI 69 [ start ]))) -1 (nil)
    (nil))

(insn 54 53 55 3 ./CppStatUtilities.cc:181 (parallel [
            (set (reg:DI 77)
                (ashift:DI (reg:DI 76)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 55 54 56 3 ./CppStatUtilities.cc:181 (parallel [
            (set (reg:DI 60 [ ivtmp.1077 ])
                (plus:DI (reg:DI 77)
                    (reg/v/f:DI 68 [ a ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 56 55 57 3 ./CppStatUtilities.cc:181 (set (reg:DI 61 [ ivtmp.1074 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 57 56 58 4 993 "" [1 uses])

(note 58 57 59 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 59 58 60 4 ("./CppStatUtilities.cc") 182)

(insn 60 59 61 4 ./CppStatUtilities.cc:182 (set (reg:DF 79)
        (mem:DF (reg:DI 60 [ ivtmp.1077 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 61 60 62 4 ./CppStatUtilities.cc:182 (set (reg:DF 78)
        (minus:DF (reg:DF 79)
            (reg/v:DF 63 [ med ]))) -1 (nil)
    (nil))

(insn 62 61 63 4 ./CppStatUtilities.cc:182 (set (reg:V2DF 80)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC33") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 63 62 64 4 ./CppStatUtilities.cc:182 (parallel [
            (set (reg:DF 81)
                (abs:DF (reg:DF 78)))
            (use (reg:V2DF 80))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 64 63 65 4 ./CppStatUtilities.cc:182 (set (mem:DF (plus:DI (mult:DI (reg:DI 61 [ ivtmp.1074 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 64 [ res ])) [3 S8 A64])
        (reg:DF 81)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 78))
        (nil)))

(insn 65 64 66 4 ./CppStatUtilities.cc:182 (parallel [
            (set (reg:DI 61 [ ivtmp.1074 ])
                (plus:DI (reg:DI 61 [ ivtmp.1074 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 66 65 67 4 ./CppStatUtilities.cc:182 (parallel [
            (set (reg:DI 60 [ ivtmp.1077 ])
                (plus:DI (reg:DI 60 [ ivtmp.1077 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(note 67 66 68 4 ("./CppStatUtilities.cc") 181)

(insn 68 67 69 4 ./CppStatUtilities.cc:181 (parallel [
            (set (reg:SI 82)
                (plus:SI (reg/v:SI 69 [ start ])
                    (subreg:SI (reg:DI 61 [ ivtmp.1074 ]) 0)))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 69 68 70 4 ./CppStatUtilities.cc:181 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 59 [ end.1080 ])
            (reg:SI 82))) -1 (nil)
    (nil))

(jump_insn 70 69 71 4 ./CppStatUtilities.cc:181 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 57)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 71 70 72 5 991 "" [1 uses])

(note 72 71 73 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(note 73 72 74 5 ("./CppStatUtilities.cc") 184)

(insn 74 73 75 5 ./CppStatUtilities.cc:184 (set (reg:SI 4 si)
        (reg/v:SI 65 [ len ])) -1 (nil)
    (nil))

(insn 75 74 76 5 ./CppStatUtilities.cc:184 (set (reg:DI 5 di)
        (reg/v/f:DI 64 [ res ])) -1 (nil)
    (nil))

(call_insn 76 75 77 5 ./CppStatUtilities.cc:184 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 77 76 78 5 ./CppStatUtilities.cc:184 (set (reg/v:DF 62 [ mad ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 78 77 79 5 ("./CppStatUtilities.cc") 186)

(insn 79 78 80 5 ./CppStatUtilities.cc:186 (set (reg:DI 5 di)
        (reg/v/f:DI 64 [ res ])) -1 (nil)
    (nil))

(call_insn 80 79 81 5 ./CppStatUtilities.cc:186 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 81 80 82 5 ./CppStatUtilities.cc:186 (set (reg:DF 83)
        (mult:DF (reg/v:DF 62 [ mad ])
            (reg/v:DF 71 [ constant ]))) -1 (nil)
    (nil))

(insn 82 81 85 5 ./CppStatUtilities.cc:186 (set (reg:DF 67 [ <result> ])
        (reg:DF 83)) -1 (nil)
    (nil))

(note 85 82 86 5 NOTE_INSN_FUNCTION_END)

(note 86 85 88 5 ("./CppStatUtilities.cc") 189)

(insn 88 86 94 5 ./CppStatUtilities.cc:189 (set (reg/i:DF 21 xmm0)
        (reg:DF 67 [ <result> ])) -1 (nil)
    (nil))

(insn 94 88 0 5 ./CppStatUtilities.cc:189 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)


;; Function void summaryStats(double*, int, summaryStats_t&, int, double) (_Z12summaryStatsPdiR14summaryStats_tid)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 3.
Deleted label in block 5.
Deleted label in block 6.
Deleted label in block 12.
Deleted label in block 14.
Forwarding edge 14->15 to 28 failed.
Deleted label in block 19.
Deleted label in block 23.
Redirecting jump 280 from 29 to 30.
Deleting block 29.


try_optimize_cfg iteration 2

Forwarding edge 14->15 to 28 failed.


try_optimize_cfg iteration 1

Forwarding edge 13->14 to 27 failed.
(note 1 0 10 ("./CppStatUtilities.cc") 69)

;; Start of basic block 0, registers live: (nil)
(note 10 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 10 4 0 ./CppStatUtilities.cc:69 (set (reg/v/f:DI 82 [ data ])
        (reg:DI 5 di [ data ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./CppStatUtilities.cc:69 (set (reg/v:SI 83 [ dataSize ])
        (reg:SI 4 si [ dataSize ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./CppStatUtilities.cc:69 (set (reg/v/f:DI 84 [ summary ])
        (reg:DI 1 dx [ summary ])) -1 (nil)
    (nil))

(insn 6 5 7 0 ./CppStatUtilities.cc:69 (set (reg/v:SI 85 [ sample ])
        (reg:SI 2 cx [ sample ])) -1 (nil)
    (nil))

(insn 7 6 8 0 ./CppStatUtilities.cc:69 (set (reg/v:DF 86 [ constant ])
        (reg:DF 21 xmm0 [ constant ])) -1 (nil)
    (nil))

(note 8 7 12 0 NOTE_INSN_FUNCTION_BEG)

(note 12 8 13 0 ("./CppStatUtilities.cc") 72)

(insn 13 12 14 0 ./CppStatUtilities.cc:72 (set (reg:CCGC 17 flags)
        (compare:CCGC (reg/v:SI 85 [ sample ])
            (reg/v:SI 83 [ dataSize ]))) -1 (nil)
    (nil))

(jump_insn 14 13 16 0 ./CppStatUtilities.cc:72 (set (pc)
        (if_then_else (gt (reg:CCGC 17 flags)
                (const_int 0 [0x0]))
            (label_ref 282)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5000 [0x1388])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 16 14 17 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(note 17 16 18 1 ("./CppStatUtilities.cc") 74)

(insn 18 17 19 1 ./CppStatUtilities.cc:74 (set (reg:DF 87)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC34") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 19 18 20 1 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -48 [0xffffffffffffffd0])) [3 prob+0 S8 A64])
        (reg:DF 87)) -1 (nil)
    (nil))

(insn 20 19 21 1 ./CppStatUtilities.cc:74 (set (reg:DF 88)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC35") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 21 20 22 1 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -40 [0xffffffffffffffd8])) [3 prob+8 S8 A64])
        (reg:DF 88)) -1 (nil)
    (nil))

(insn 22 21 23 1 ./CppStatUtilities.cc:74 (set (reg:DF 89)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC36") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 23 22 24 1 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -32 [0xffffffffffffffe0])) [3 prob+16 S8 A64])
        (reg:DF 89)) -1 (nil)
    (nil))

(insn 24 23 25 1 ./CppStatUtilities.cc:74 (set (reg:DF 90)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC37") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 25 24 26 1 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -24 [0xffffffffffffffe8])) [3 prob+24 S8 A64])
        (reg:DF 90)) -1 (nil)
    (nil))

(insn 26 25 27 1 ./CppStatUtilities.cc:74 (set (reg:DF 91)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC38") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 27 26 28 1 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -16 [0xfffffffffffffff0])) [3 prob+32 S8 A64])
        (reg:DF 91)) -1 (nil)
    (nil))

(note 28 27 29 1 ("./CppStatUtilities.cc") 77)

(insn 29 28 30 1 ./CppStatUtilities.cc:77 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 85 [ sample ])
            (reg/v:SI 83 [ dataSize ]))) -1 (nil)
    (nil))

(jump_insn 30 29 32 1 ./CppStatUtilities.cc:77 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 112)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 5120 [0x1400])
        (nil)))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(note 32 30 33 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 33 32 34 2 ./CppStatUtilities.cc:77 (set (reg/v:SI 75 [ xLen ])
        (reg/v:SI 85 [ sample ])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 34 33 35 3 1034 "" [1 uses])

(note 35 34 36 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 36 35 37 3 ("./CppStatUtilities.cc") 79)

(insn 37 36 38 3 ./CppStatUtilities.cc:79 (parallel [
            (set (reg:DI 92)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 38 37 39 3 ./CppStatUtilities.cc:79 (set (reg:SI 37 r8)
        (const_int 100 [0x64])) -1 (nil)
    (nil))

(insn 39 38 40 3 ./CppStatUtilities.cc:79 (set (reg:SI 2 cx)
        (const_int 5 [0x5])) -1 (nil)
    (nil))

(insn 40 39 41 3 ./CppStatUtilities.cc:79 (set (reg:DI 1 dx)
        (reg:DI 92)) -1 (nil)
    (nil))

(insn 41 40 42 3 ./CppStatUtilities.cc:79 (set (reg:SI 4 si)
        (reg/v:SI 75 [ xLen ])) -1 (nil)
    (nil))

(insn 42 41 43 3 ./CppStatUtilities.cc:79 (set (reg:DI 5 di)
        (reg/v/f:DI 82 [ data ])) -1 (nil)
    (nil))

(call_insn 43 42 44 3 ./CppStatUtilities.cc:79 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Z11percentilesPKdiS0_ii") [flags 0x3] <function_decl 0x2b5eb636a000 percentiles>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 37 r8))
                        (nil)))))))

(insn 44 43 45 3 ./CppStatUtilities.cc:79 (set (reg/v/f:DI 74 [ perc ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(note 45 44 46 3 ("./CppStatUtilities.cc") 81)

(insn 46 45 47 3 ./CppStatUtilities.cc:81 (set (reg:DI 93)
        (sign_extend:DI (reg/v:SI 75 [ xLen ]))) -1 (nil)
    (nil))

(insn 47 46 48 3 ./CppStatUtilities.cc:81 (parallel [
            (set (reg:DI 94)
                (ashift:DI (reg:DI 93)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 48 47 49 3 ./CppStatUtilities.cc:81 (set (reg:DI 5 di)
        (reg:DI 94)) -1 (nil)
    (nil))

(call_insn 49 48 50 3 ./CppStatUtilities.cc:81 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 50 49 51 3 ./CppStatUtilities.cc:81 (set (reg/f:DI 95)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 95)
        (nil)))

(insn 51 50 52 3 ./CppStatUtilities.cc:81 (set (reg:DI 80 [ D.34502 ])
        (reg/f:DI 95)) -1 (nil)
    (nil))

(note 52 51 53 3 ("./CppStatUtilities.cc") 82)

(insn 53 52 54 3 ./CppStatUtilities.cc:82 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 80 [ D.34502 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 54 53 56 3 ./CppStatUtilities.cc:82 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 67)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 100 [0x64])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 56 54 57 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(note 57 56 58 4 ("./CppStatUtilities.cc") 81)

(insn 58 57 59 4 ./CppStatUtilities.cc:81 (set (reg/v/f:DI 72 [ res ])
        (reg:DI 80 [ D.34502 ])) -1 (nil)
    (nil))

(note 59 58 60 4 ("./CppStatUtilities.cc") 84)

(insn 60 59 61 4 ./CppStatUtilities.cc:84 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 75 [ xLen ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 61 60 63 4 ./CppStatUtilities.cc:84 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 79)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(note 63 61 64 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 64 63 65 5 ./CppStatUtilities.cc:84 (parallel [
            (set (reg:DI 65 [ prephitmp.1125 ])
                (plus:DI (reg/v/f:DI 74 [ perc ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 65 64 66 5 ./CppStatUtilities.cc:84 (set (pc)
        (label_ref 99)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
 (nil)

(barrier 66 65 67)

;; Start of basic block 6, registers live: (nil)
(code_label 67 66 68 6 1035 "" [1 uses])

(note 68 67 69 6 [bb 6] NOTE_INSN_BASIC_BLOCK)

(note 69 68 70 6 ("./CppStatUtilities.cc") 82)

(insn 70 69 71 6 ./CppStatUtilities.cc:82 (set (reg:SI 2 cx)
        (const_int 82 [0x52])) -1 (nil)
    (nil))

(insn 71 70 72 6 ./CppStatUtilities.cc:82 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 72 71 73 6 ./CppStatUtilities.cc:82 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 73 72 74 6 ./CppStatUtilities.cc:82 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 74 73 75 6 ./CppStatUtilities.cc:82 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 75 74 76 6 ./CppStatUtilities.cc:82 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 76 75 77 6 ./CppStatUtilities.cc:82 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 77 76 78 6 ./CppStatUtilities.cc:82 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 6, registers live:
 (nil)

(barrier 78 77 79)

;; Start of basic block 7, registers live: (nil)
(code_label 79 78 80 7 1037 "" [1 uses])

(note 80 79 81 7 [bb 7] NOTE_INSN_BASIC_BLOCK)

(insn 81 80 82 7 ./CppStatUtilities.cc:82 (parallel [
            (set (reg:DI 65 [ prephitmp.1125 ])
                (plus:DI (reg/v/f:DI 74 [ perc ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 82 81 83 7 ./CppStatUtilities.cc:82 (set (reg/v:SI 71 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 83 82 84 7 ./CppStatUtilities.cc:82 (set (reg:DI 63 [ ivtmp.1144 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 7, registers live:
 (nil)

;; Start of basic block 8, registers live: (nil)
(code_label 84 83 85 8 1040 "" [1 uses])

(note 85 84 86 8 [bb 8] NOTE_INSN_BASIC_BLOCK)

(note 86 85 87 8 ("./CppStatUtilities.cc") 69)

(insn 87 86 88 8 ./CppStatUtilities.cc:69 (set (reg:DI 62 [ D.37889 ])
        (reg:DI 63 [ ivtmp.1144 ])) -1 (nil)
    (nil))

(note 88 87 89 8 ("./CppStatUtilities.cc") 85)

(insn 89 88 90 8 ./CppStatUtilities.cc:85 (set (reg:DF 97)
        (mem:DF (plus:DI (mult:DI (reg:DI 62 [ D.37889 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 82 [ data ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 90 89 91 8 ./CppStatUtilities.cc:85 (set (reg:DF 96)
        (minus:DF (reg:DF 97)
            (mem:DF (reg:DI 65 [ prephitmp.1125 ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 91 90 92 8 ./CppStatUtilities.cc:85 (set (reg:V2DF 98)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC39") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 92 91 93 8 ./CppStatUtilities.cc:85 (parallel [
            (set (reg:DF 99)
                (abs:DF (reg:DF 96)))
            (use (reg:V2DF 98))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 93 92 94 8 ./CppStatUtilities.cc:85 (set (mem:DF (plus:DI (mult:DI (reg:DI 62 [ D.37889 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 72 [ res ])) [3 S8 A64])
        (reg:DF 99)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 96))
        (nil)))

(note 94 93 95 8 ("./CppStatUtilities.cc") 84)

(insn 95 94 96 8 ./CppStatUtilities.cc:84 (parallel [
            (set (reg/v:SI 71 [ i ])
                (plus:SI (reg/v:SI 71 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 96 95 97 8 ./CppStatUtilities.cc:84 (parallel [
            (set (reg:DI 63 [ ivtmp.1144 ])
                (plus:DI (reg:DI 63 [ ivtmp.1144 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 97 96 98 8 ./CppStatUtilities.cc:84 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 71 [ i ])
            (reg/v:SI 75 [ xLen ]))) -1 (nil)
    (nil))

(jump_insn 98 97 99 8 ./CppStatUtilities.cc:84 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 84)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 8, registers live:
 (nil)

;; Start of basic block 9, registers live: (nil)
(code_label 99 98 100 9 1039 "" [1 uses])

(note 100 99 101 9 [bb 9] NOTE_INSN_BASIC_BLOCK)

(note 101 100 102 9 ("./CppStatUtilities.cc") 87)

(insn 102 101 103 9 ./CppStatUtilities.cc:87 (set (reg:SI 4 si)
        (reg/v:SI 75 [ xLen ])) -1 (nil)
    (nil))

(insn 103 102 104 9 ./CppStatUtilities.cc:87 (set (reg:DI 5 di)
        (reg/v/f:DI 72 [ res ])) -1 (nil)
    (nil))

(call_insn 104 103 105 9 ./CppStatUtilities.cc:87 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 105 104 106 9 ./CppStatUtilities.cc:87 (set (reg:DF 79 [ D.34520 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 106 105 107 9 ./CppStatUtilities.cc:87 (set (reg/v:DF 73 [ mad ])
        (mult:DF (reg:DF 79 [ D.34520 ])
            (reg/v:DF 86 [ constant ]))) -1 (nil)
    (nil))

(note 107 106 108 9 ("./CppStatUtilities.cc") 89)

(insn 108 107 109 9 ./CppStatUtilities.cc:89 (set (reg:DI 5 di)
        (reg/v/f:DI 72 [ res ])) -1 (nil)
    (nil))

(call_insn 109 108 110 9 ./CppStatUtilities.cc:89 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(jump_insn 110 109 111 9 ./CppStatUtilities.cc:89 (set (pc)
        (label_ref 229)) -1 (nil)
    (nil))
;; End of basic block 9, registers live:
 (nil)

(barrier 111 110 112)

;; Start of basic block 10, registers live: (nil)
(code_label 112 111 113 10 1032 "" [1 uses])

(note 113 112 114 10 [bb 10] NOTE_INSN_BASIC_BLOCK)

(note 114 113 115 10 ("./CppStatUtilities.cc") 93)

(insn 115 114 116 10 ./CppStatUtilities.cc:93 (set (reg:DI 100)
        (sign_extend:DI (reg/v:SI 85 [ sample ]))) -1 (nil)
    (nil))

(insn 116 115 117 10 ./CppStatUtilities.cc:93 (parallel [
            (set (reg:DI 81 [ D.34501 ])
                (ashift:DI (reg:DI 100)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 117 116 118 10 ./CppStatUtilities.cc:93 (set (reg:DI 5 di)
        (reg:DI 81 [ D.34501 ])) -1 (nil)
    (nil))

(call_insn 118 117 119 10 ./CppStatUtilities.cc:93 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 119 118 120 10 ./CppStatUtilities.cc:93 (set (reg/f:DI 101)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 101)
        (nil)))

(insn 120 119 121 10 ./CppStatUtilities.cc:93 (set (reg:DI 59 [ ivtmp.1158 ])
        (reg/f:DI 101)) -1 (nil)
    (nil))

(note 121 120 122 10 ("./CppStatUtilities.cc") 94)

(insn 122 121 123 10 ./CppStatUtilities.cc:94 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 59 [ ivtmp.1158 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 123 122 125 10 ./CppStatUtilities.cc:94 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 135)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9900 [0x26ac])
        (nil)))
;; End of basic block 10, registers live:
 (nil)

;; Start of basic block 11, registers live: (nil)
(note 125 123 126 11 [bb 11] NOTE_INSN_BASIC_BLOCK)

(insn 126 125 127 11 ./CppStatUtilities.cc:94 (set (reg:SI 2 cx)
        (const_int 94 [0x5e])) -1 (nil)
    (nil))

(insn 127 126 128 11 ./CppStatUtilities.cc:94 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 128 127 129 11 ./CppStatUtilities.cc:94 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 129 128 130 11 ./CppStatUtilities.cc:94 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 130 129 131 11 ./CppStatUtilities.cc:94 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 131 130 132 11 ./CppStatUtilities.cc:94 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 132 131 133 11 ./CppStatUtilities.cc:94 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 133 132 134 11 ./CppStatUtilities.cc:94 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 11, registers live:
 (nil)

(barrier 134 133 135)

;; Start of basic block 12, registers live: (nil)
(code_label 135 134 136 12 1042 "" [1 uses])

(note 136 135 137 12 [bb 12] NOTE_INSN_BASIC_BLOCK)

(note 137 136 138 12 ("./CppStatUtilities.cc") 96)

(insn 138 137 139 12 ./CppStatUtilities.cc:96 (set (reg:DI 5 di)
        (reg:DI 81 [ D.34501 ])) -1 (nil)
    (nil))

(call_insn 139 138 140 12 ./CppStatUtilities.cc:96 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b5eb3cce100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 140 139 141 12 ./CppStatUtilities.cc:96 (set (reg/f:DI 102)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 102)
        (nil)))

(insn 141 140 142 12 ./CppStatUtilities.cc:96 (set (reg:DI 78 [ D.34524 ])
        (reg/f:DI 102)) -1 (nil)
    (nil))

(note 142 141 143 12 ("./CppStatUtilities.cc") 97)

(insn 143 142 144 12 ./CppStatUtilities.cc:97 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg:DI 78 [ D.34524 ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 144 143 146 12 ./CppStatUtilities.cc:97 (set (pc)
        (if_then_else (eq (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 157)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 100 [0x64])
        (nil)))
;; End of basic block 12, registers live:
 (nil)

;; Start of basic block 13, registers live: (nil)
(note 146 144 147 13 [bb 13] NOTE_INSN_BASIC_BLOCK)

(note 147 146 148 13 ("./CppStatUtilities.cc") 93)

(insn 148 147 149 13 ./CppStatUtilities.cc:93 (set (reg/v/f:DI 70 [ x ])
        (reg:DI 59 [ ivtmp.1158 ])) -1 (nil)
    (nil))

(note 149 148 150 13 ("./CppStatUtilities.cc") 96)

(insn 150 149 151 13 ./CppStatUtilities.cc:96 (set (reg/v/f:DI 69 [ res ])
        (reg:DI 78 [ D.34524 ])) -1 (nil)
    (nil))

(note 151 150 152 13 ("./CppStatUtilities.cc") 100)

(insn 152 151 153 13 ./CppStatUtilities.cc:100 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 85 [ sample ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 153 152 156 13 ./CppStatUtilities.cc:100 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 169)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 13, registers live:
 (nil)

;; Start of basic block 14, registers live: (nil)
(note 156 153 154 14 [bb 14] NOTE_INSN_BASIC_BLOCK)

(jump_insn 154 156 155 14 ./CppStatUtilities.cc:100 (set (pc)
        (label_ref 298)) -1 (nil)
    (nil))
;; End of basic block 14, registers live:
 (nil)

(barrier 155 154 157)

;; Start of basic block 15, registers live: (nil)
(code_label 157 155 158 15 1044 "" [1 uses])

(note 158 157 159 15 [bb 15] NOTE_INSN_BASIC_BLOCK)

(note 159 158 160 15 ("./CppStatUtilities.cc") 97)

(insn 160 159 161 15 ./CppStatUtilities.cc:97 (set (reg:SI 2 cx)
        (const_int 97 [0x61])) -1 (nil)
    (nil))

(insn 161 160 162 15 ./CppStatUtilities.cc:97 (set (reg:DI 1 dx)
        (symbol_ref/f:DI ("*.LC3") [flags 0x2] <string_cst 0x2b5eb671c300>)) -1 (nil)
    (nil))

(insn 162 161 163 15 ./CppStatUtilities.cc:97 (set (reg:DI 4 si)
        (symbol_ref/f:DI ("*.LC4") [flags 0x2] <string_cst 0x2b5eb66c6f50>)) -1 (nil)
    (nil))

(insn 163 162 164 15 ./CppStatUtilities.cc:97 (set (reg:DI 5 di)
        (mem/f/c/i:DI (symbol_ref:DI ("stderr") [flags 0x40] <var_decl 0x2b5eb42da790 stderr>) [11 stderr+0 S8 A64])) -1 (nil)
    (nil))

(insn 164 163 165 15 ./CppStatUtilities.cc:97 (set (reg:QI 0 ax)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(call_insn 165 164 166 15 ./CppStatUtilities.cc:97 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("fprintf") [flags 0x41] <function_decl 0x2b5eb3c9fd00 fprintf>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:QI 0 ax))
        (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
            (expr_list:REG_DEP_TRUE (use (reg:DI 4 si))
                (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                        (nil)))))))

(insn 166 165 167 15 ./CppStatUtilities.cc:97 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn 167 166 168 15 ./CppStatUtilities.cc:97 (call (mem:QI (symbol_ref:DI ("exit") [flags 0x41] <function_decl 0x2b5eb3cc3b00 exit>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_NORETURN (const_int 0 [0x0])
        (expr_list:REG_EH_REGION (const_int 0 [0x0])
            (nil)))
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (nil)))
;; End of basic block 15, registers live:
 (nil)

(barrier 168 167 169)

;; Start of basic block 16, registers live: (nil)
(code_label 169 168 170 16 1046 "" [1 uses])

(note 170 169 171 16 [bb 16] NOTE_INSN_BASIC_BLOCK)

(insn 171 170 172 16 ./CppStatUtilities.cc:97 (set (reg/v:SI 68 [ k ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 16, registers live:
 (nil)

;; Start of basic block 17, registers live: (nil)
(code_label 172 171 173 17 1048 "" [1 uses])

(note 173 172 174 17 [bb 17] NOTE_INSN_BASIC_BLOCK)

(note 174 173 175 17 ("./CppStatUtilities.cc") 101)

(call_insn 175 174 176 17 ./CppStatUtilities.cc:101 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("rand") [flags 0x41] <function_decl 0x2b5eb50eef00 rand>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (nil))

(insn 176 175 177 17 ./CppStatUtilities.cc:101 (set (reg:SI 77 [ D.34535 ])
        (reg:SI 0 ax)) -1 (nil)
    (nil))

(insn 177 176 178 17 ./CppStatUtilities.cc:101 (parallel [
            (set (reg:SI 105)
                (div:SI (reg:SI 77 [ D.34535 ])
                    (reg/v:SI 85 [ sample ])))
            (set (reg:SI 104)
                (mod:SI (reg:SI 77 [ D.34535 ])
                    (reg/v:SI 85 [ sample ])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 178 177 179 17 ./CppStatUtilities.cc:101 (set (reg:DI 106)
        (sign_extend:DI (reg:SI 104))) -1 (nil)
    (nil))

(insn 179 178 180 17 ./CppStatUtilities.cc:101 (set (reg:DF 107)
        (mem:DF (plus:DI (mult:DI (reg:DI 106)
                    (const_int 8 [0x8]))
                (reg/v/f:DI 82 [ data ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 180 179 181 17 ./CppStatUtilities.cc:101 (set (mem:DF (reg:DI 59 [ ivtmp.1158 ]) [3 S8 A64])
        (reg:DF 107)) -1 (nil)
    (nil))

(note 181 180 182 17 ("./CppStatUtilities.cc") 100)

(insn 182 181 183 17 ./CppStatUtilities.cc:100 (parallel [
            (set (reg/v:SI 68 [ k ])
                (plus:SI (reg/v:SI 68 [ k ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 183 182 184 17 ./CppStatUtilities.cc:100 (parallel [
            (set (reg:DI 59 [ ivtmp.1158 ])
                (plus:DI (reg:DI 59 [ ivtmp.1158 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 184 183 185 17 ./CppStatUtilities.cc:100 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 68 [ k ])
            (reg/v:SI 85 [ sample ]))) -1 (nil)
    (nil))

(jump_insn 185 184 187 17 ./CppStatUtilities.cc:100 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 172)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 17, registers live:
 (nil)

;; Start of basic block 18, registers live: (nil)
(note 187 185 188 18 [bb 18] NOTE_INSN_BASIC_BLOCK)

(note 188 187 189 18 ("./CppStatUtilities.cc") 103)

(insn 189 188 190 18 ./CppStatUtilities.cc:103 (parallel [
            (set (reg:DI 108)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 190 189 191 18 ./CppStatUtilities.cc:103 (set (reg:SI 37 r8)
        (const_int 100 [0x64])) -1 (nil)
    (nil))

(insn 191 190 192 18 ./CppStatUtilities.cc:103 (set (reg:SI 2 cx)
        (const_int 5 [0x5])) -1 (nil)
    (nil))

(insn 192 191 193 18 ./CppStatUtilities.cc:103 (set (reg:DI 1 dx)
        (reg:DI 108)) -1 (nil)
    (nil))

(insn 193 192 194 18 ./CppStatUtilities.cc:103 (set (reg:SI 4 si)
        (reg/v:SI 85 [ sample ])) -1 (nil)
    (nil))

(insn 194 193 195 18 ./CppStatUtilities.cc:103 (set (reg:DI 5 di)
        (reg/v/f:DI 70 [ x ])) -1 (nil)
    (nil))

(call_insn 195 194 196 18 ./CppStatUtilities.cc:103 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Z11percentilesPKdiS0_ii") [flags 0x3] <function_decl 0x2b5eb636a000 percentiles>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 37 r8))
                        (nil)))))))

(insn 196 195 197 18 ./CppStatUtilities.cc:103 (set (reg/v/f:DI 74 [ perc ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 197 196 198 18 ./CppStatUtilities.cc:103 (parallel [
            (set (reg:DI 65 [ prephitmp.1125 ])
                (plus:DI (reg/v/f:DI 74 [ perc ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 198 197 199 18 ./CppStatUtilities.cc:103 (set (reg/v:SI 67 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 199 198 200 18 ./CppStatUtilities.cc:103 (set (reg:DI 61 [ ivtmp.1150 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 18, registers live:
 (nil)

;; Start of basic block 19, registers live: (nil)
(code_label 200 199 201 19 1050 "" [1 uses])

(note 201 200 202 19 [bb 19] NOTE_INSN_BASIC_BLOCK)

(note 202 201 203 19 ("./CppStatUtilities.cc") 69)

(insn 203 202 204 19 ./CppStatUtilities.cc:69 (set (reg:DI 60 [ D.37899 ])
        (reg:DI 61 [ ivtmp.1150 ])) -1 (nil)
    (nil))

(note 204 203 205 19 ("./CppStatUtilities.cc") 106)

(insn 205 204 206 19 ./CppStatUtilities.cc:106 (set (reg:DF 110)
        (mem:DF (plus:DI (mult:DI (reg:DI 60 [ D.37899 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 70 [ x ])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 206 205 207 19 ./CppStatUtilities.cc:106 (set (reg:DF 109)
        (minus:DF (reg:DF 110)
            (mem:DF (reg:DI 65 [ prephitmp.1125 ]) [3 S8 A64]))) -1 (nil)
    (nil))

(insn 207 206 208 19 ./CppStatUtilities.cc:106 (set (reg:V2DF 111)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC39") [flags 0x2]) [3 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 208 207 209 19 ./CppStatUtilities.cc:106 (parallel [
            (set (reg:DF 112)
                (abs:DF (reg:DF 109)))
            (use (reg:V2DF 111))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 209 208 210 19 ./CppStatUtilities.cc:106 (set (mem:DF (plus:DI (mult:DI (reg:DI 60 [ D.37899 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 69 [ res ])) [3 S8 A64])
        (reg:DF 112)) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 109))
        (nil)))

(note 210 209 211 19 ("./CppStatUtilities.cc") 105)

(insn 211 210 212 19 ./CppStatUtilities.cc:105 (parallel [
            (set (reg/v:SI 67 [ i ])
                (plus:SI (reg/v:SI 67 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 212 211 213 19 ./CppStatUtilities.cc:105 (parallel [
            (set (reg:DI 61 [ ivtmp.1150 ])
                (plus:DI (reg:DI 61 [ ivtmp.1150 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 213 212 214 19 ./CppStatUtilities.cc:105 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 67 [ i ])
            (reg/v:SI 85 [ sample ]))) -1 (nil)
    (nil))

(jump_insn 214 213 215 19 ./CppStatUtilities.cc:105 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 200)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 19, registers live:
 (nil)

;; Start of basic block 20, registers live: (nil)
(code_label 215 214 216 20 1051 "" [1 uses])

(note 216 215 217 20 [bb 20] NOTE_INSN_BASIC_BLOCK)

(note 217 216 218 20 ("./CppStatUtilities.cc") 108)

(insn 218 217 219 20 ./CppStatUtilities.cc:108 (set (reg:SI 4 si)
        (reg/v:SI 85 [ sample ])) -1 (nil)
    (nil))

(insn 219 218 220 20 ./CppStatUtilities.cc:108 (set (reg:DI 5 di)
        (reg/v/f:DI 69 [ res ])) -1 (nil)
    (nil))

(call_insn 220 219 221 20 ./CppStatUtilities.cc:108 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("medianInPlace") [flags 0x41] <function_decl 0x2b5eb631e600 medianInPlace>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 221 220 222 20 ./CppStatUtilities.cc:108 (set (reg:DF 76 [ D.34555 ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(insn 222 221 223 20 ./CppStatUtilities.cc:108 (set (reg/v:DF 73 [ mad ])
        (mult:DF (reg:DF 76 [ D.34555 ])
            (reg/v:DF 86 [ constant ]))) -1 (nil)
    (nil))

(note 223 222 224 20 ("./CppStatUtilities.cc") 110)

(insn 224 223 225 20 ./CppStatUtilities.cc:110 (set (reg:DI 5 di)
        (reg/v/f:DI 69 [ res ])) -1 (nil)
    (nil))

(call_insn 225 224 226 20 ./CppStatUtilities.cc:110 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 226 225 227 20 ("./CppStatUtilities.cc") 111)

(insn 227 226 228 20 ./CppStatUtilities.cc:111 (set (reg:DI 5 di)
        (reg/v/f:DI 70 [ x ])) -1 (nil)
    (nil))

(call_insn 228 227 229 20 ./CppStatUtilities.cc:111 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 20, registers live:
 (nil)

;; Start of basic block 21, registers live: (nil)
(code_label 229 228 230 21 1041 "" [1 uses])

(note 230 229 231 21 [bb 21] NOTE_INSN_BASIC_BLOCK)

(note 231 230 232 21 ("./CppStatUtilities.cc") 114)

(insn 232 231 233 21 ./CppStatUtilities.cc:114 (set (reg:DF 113)
        (mem:DF (reg/v/f:DI 74 [ perc ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 233 232 234 21 ./CppStatUtilities.cc:114 (set (mem/s:DF (reg/v/f:DI 84 [ summary ]) [3 <variable>.min+0 S8 A64])
        (reg:DF 113)) -1 (nil)
    (nil))

(note 234 233 235 21 ("./CppStatUtilities.cc") 115)

(insn 235 234 236 21 ./CppStatUtilities.cc:115 (set (reg:DF 114)
        (mem:DF (plus:DI (reg/v/f:DI 74 [ perc ])
                (const_int 8 [0x8])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 236 235 237 21 ./CppStatUtilities.cc:115 (set (mem/s:DF (plus:DI (reg/v/f:DI 84 [ summary ])
                (const_int 8 [0x8])) [3 <variable>.max+0 S8 A64])
        (reg:DF 114)) -1 (nil)
    (nil))

(note 237 236 238 21 ("./CppStatUtilities.cc") 116)

(insn 238 237 239 21 ./CppStatUtilities.cc:116 (set (reg:DF 115)
        (mem:DF (reg:DI 65 [ prephitmp.1125 ]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 239 238 240 21 ./CppStatUtilities.cc:116 (set (mem/s:DF (plus:DI (reg/v/f:DI 84 [ summary ])
                (const_int 16 [0x10])) [3 <variable>.med+0 S8 A64])
        (reg:DF 115)) -1 (nil)
    (nil))

(note 240 239 241 21 ("./CppStatUtilities.cc") 117)

(insn 241 240 242 21 ./CppStatUtilities.cc:117 (set (mem/s:DF (plus:DI (reg/v/f:DI 84 [ summary ])
                (const_int 32 [0x20])) [3 <variable>.mad+0 S8 A64])
        (reg/v:DF 73 [ mad ])) -1 (nil)
    (nil))

(note 242 241 243 21 ("./CppStatUtilities.cc") 118)

(insn 243 242 244 21 ./CppStatUtilities.cc:118 (set (reg:DF 116)
        (mem:DF (plus:DI (reg/v/f:DI 74 [ perc ])
                (const_int 24 [0x18])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 244 243 245 21 ./CppStatUtilities.cc:118 (set (mem/s:DF (plus:DI (reg/v/f:DI 84 [ summary ])
                (const_int 40 [0x28])) [3 <variable>.q1+0 S8 A64])
        (reg:DF 116)) -1 (nil)
    (nil))

(note 245 244 246 21 ("./CppStatUtilities.cc") 119)

(insn 246 245 247 21 ./CppStatUtilities.cc:119 (set (reg:DF 117)
        (mem:DF (plus:DI (reg/v/f:DI 74 [ perc ])
                (const_int 32 [0x20])) [3 S8 A64])) -1 (nil)
    (nil))

(insn 247 246 248 21 ./CppStatUtilities.cc:119 (set (mem/s:DF (plus:DI (reg/v/f:DI 84 [ summary ])
                (const_int 48 [0x30])) [3 <variable>.q3+0 S8 A64])
        (reg:DF 117)) -1 (nil)
    (nil))

(note 248 247 249 21 ("./CppStatUtilities.cc") 123)

(insn 249 248 250 21 ./CppStatUtilities.cc:123 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 83 [ dataSize ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 250 249 252 21 ./CppStatUtilities.cc:123 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 256)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 21, registers live:
 (nil)

;; Start of basic block 22, registers live: (nil)
(note 252 250 253 22 [bb 22] NOTE_INSN_BASIC_BLOCK)

(insn 253 252 254 22 ./CppStatUtilities.cc:123 (set (reg/v:DF 58 [ mean.1166 ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC34") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(jump_insn 254 253 255 22 ./CppStatUtilities.cc:123 (set (pc)
        (label_ref 270)) -1 (nil)
    (nil))
;; End of basic block 22, registers live:
 (nil)

(barrier 255 254 256)

;; Start of basic block 23, registers live: (nil)
(code_label 256 255 257 23 1052 "" [1 uses])

(note 257 256 258 23 [bb 23] NOTE_INSN_BASIC_BLOCK)

(insn 258 257 259 23 ./CppStatUtilities.cc:123 (set (reg:DI 64 [ ivtmp.1141 ])
        (reg/v/f:DI 82 [ data ])) -1 (nil)
    (nil))

(insn 259 258 260 23 ./CppStatUtilities.cc:123 (set (reg/v:DF 58 [ mean.1166 ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC34") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(insn 260 259 261 23 ./CppStatUtilities.cc:123 (set (reg/v:SI 66 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 23, registers live:
 (nil)

;; Start of basic block 24, registers live: (nil)
(code_label 261 260 262 24 1055 "" [1 uses])

(note 262 261 263 24 [bb 24] NOTE_INSN_BASIC_BLOCK)

(note 263 262 264 24 ("./CppStatUtilities.cc") 124)

(insn 264 263 265 24 ./CppStatUtilities.cc:124 (set (reg/v:DF 58 [ mean.1166 ])
        (plus:DF (reg/v:DF 58 [ mean.1166 ])
            (mem:DF (reg:DI 64 [ ivtmp.1141 ]) [3 S8 A64]))) -1 (nil)
    (nil))

(note 265 264 266 24 ("./CppStatUtilities.cc") 123)

(insn 266 265 267 24 ./CppStatUtilities.cc:123 (parallel [
            (set (reg/v:SI 66 [ i ])
                (plus:SI (reg/v:SI 66 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 267 266 268 24 ./CppStatUtilities.cc:123 (parallel [
            (set (reg:DI 64 [ ivtmp.1141 ])
                (plus:DI (reg:DI 64 [ ivtmp.1141 ])
                    (const_int 8 [0x8])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 268 267 269 24 ./CppStatUtilities.cc:123 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 66 [ i ])
            (reg/v:SI 83 [ dataSize ]))) -1 (nil)
    (nil))

(jump_insn 269 268 270 24 ./CppStatUtilities.cc:123 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 261)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 24, registers live:
 (nil)

;; Start of basic block 25, registers live: (nil)
(code_label 270 269 271 25 1054 "" [1 uses])

(note 271 270 272 25 [bb 25] NOTE_INSN_BASIC_BLOCK)

(note 272 271 273 25 ("./CppStatUtilities.cc") 127)

(insn 273 272 274 25 ./CppStatUtilities.cc:127 (set (reg:DF 118)
        (float:DF (reg/v:SI 83 [ dataSize ]))) -1 (nil)
    (nil))

(insn 274 273 275 25 ./CppStatUtilities.cc:127 (set (reg:DF 119)
        (div:DF (reg/v:DF 58 [ mean.1166 ])
            (reg:DF 118))) -1 (nil)
    (nil))

(insn 275 274 276 25 ./CppStatUtilities.cc:127 (set (mem/s:DF (plus:DI (reg/v/f:DI 84 [ summary ])
                (const_int 24 [0x18])) [3 <variable>.mean+0 S8 A64])
        (reg:DF 119)) -1 (nil)
    (nil))

(note 276 275 277 25 ("./CppStatUtilities.cc") 129)

(insn 277 276 278 25 ./CppStatUtilities.cc:129 (set (reg:DI 5 di)
        (reg/v/f:DI 74 [ perc ])) -1 (nil)
    (nil))

(call_insn 278 277 279 25 ./CppStatUtilities.cc:129 (call (mem:QI (symbol_ref:DI ("free") [flags 0x41] <function_decl 0x2b5eb5102400 free>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(note 279 278 280 25 ("./CppStatUtilities.cc") 130)

(jump_insn 280 279 281 25 ./CppStatUtilities.cc:130 (set (pc)
        (label_ref:DI 317)) 548 {jump} (nil)
    (nil))
;; End of basic block 25, registers live:
 (nil)

(barrier 281 280 282)

;; Start of basic block 26, registers live: (nil)
(code_label 282 281 283 26 1030 "" [1 uses])

(note 283 282 284 26 [bb 26] NOTE_INSN_BASIC_BLOCK)

(note 284 283 285 26 ("./CppStatUtilities.cc") 74)

(insn 285 284 286 26 ./CppStatUtilities.cc:74 (set (reg:DF 120)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC34") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 286 285 287 26 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -48 [0xffffffffffffffd0])) [3 prob+0 S8 A64])
        (reg:DF 120)) -1 (nil)
    (nil))

(insn 287 286 288 26 ./CppStatUtilities.cc:74 (set (reg:DF 121)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC35") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 288 287 289 26 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -40 [0xffffffffffffffd8])) [3 prob+8 S8 A64])
        (reg:DF 121)) -1 (nil)
    (nil))

(insn 289 288 290 26 ./CppStatUtilities.cc:74 (set (reg:DF 122)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC36") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 290 289 291 26 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -32 [0xffffffffffffffe0])) [3 prob+16 S8 A64])
        (reg:DF 122)) -1 (nil)
    (nil))

(insn 291 290 292 26 ./CppStatUtilities.cc:74 (set (reg:DF 123)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC37") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 292 291 293 26 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -24 [0xffffffffffffffe8])) [3 prob+24 S8 A64])
        (reg:DF 123)) -1 (nil)
    (nil))

(insn 293 292 294 26 ./CppStatUtilities.cc:74 (set (reg:DF 124)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC38") [flags 0x2]) [3 S8 A64])) -1 (nil)
    (nil))

(insn 294 293 295 26 ./CppStatUtilities.cc:74 (set (mem/s:DF (plus:DI (reg/f:DI 54 virtual-stack-vars)
                (const_int -16 [0xfffffffffffffff0])) [3 prob+32 S8 A64])
        (reg:DF 124)) -1 (nil)
    (nil))

(insn 295 294 296 26 ./CppStatUtilities.cc:74 (set (reg/v:SI 75 [ xLen ])
        (reg/v:SI 83 [ dataSize ])) -1 (nil)
    (nil))

(jump_insn 296 295 297 26 ./CppStatUtilities.cc:74 (set (pc)
        (label_ref 34)) -1 (nil)
    (nil))
;; End of basic block 26, registers live:
 (nil)

(barrier 297 296 298)

;; Start of basic block 27, registers live: (nil)
(code_label 298 297 299 27 1047 "" [1 uses])

(note 299 298 300 27 [bb 27] NOTE_INSN_BASIC_BLOCK)

(note 300 299 301 27 ("./CppStatUtilities.cc") 103)

(insn 301 300 302 27 ./CppStatUtilities.cc:103 (parallel [
            (set (reg:DI 125)
                (plus:DI (reg/f:DI 54 virtual-stack-vars)
                    (const_int -48 [0xffffffffffffffd0])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 302 301 303 27 ./CppStatUtilities.cc:103 (set (reg:SI 37 r8)
        (const_int 100 [0x64])) -1 (nil)
    (nil))

(insn 303 302 304 27 ./CppStatUtilities.cc:103 (set (reg:SI 2 cx)
        (const_int 5 [0x5])) -1 (nil)
    (nil))

(insn 304 303 305 27 ./CppStatUtilities.cc:103 (set (reg:DI 1 dx)
        (reg:DI 125)) -1 (nil)
    (nil))

(insn 305 304 306 27 ./CppStatUtilities.cc:103 (set (reg:SI 4 si)
        (reg/v:SI 85 [ sample ])) -1 (nil)
    (nil))

(insn 306 305 307 27 ./CppStatUtilities.cc:103 (set (reg:DI 5 di)
        (reg/v/f:DI 70 [ x ])) -1 (nil)
    (nil))

(call_insn 307 306 308 27 ./CppStatUtilities.cc:103 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("_Z11percentilesPKdiS0_ii") [flags 0x3] <function_decl 0x2b5eb636a000 percentiles>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (expr_list:REG_DEP_TRUE (use (reg:DI 1 dx))
                (expr_list:REG_DEP_TRUE (use (reg:SI 2 cx))
                    (expr_list:REG_DEP_TRUE (use (reg:SI 37 r8))
                        (nil)))))))

(insn 308 307 309 27 ./CppStatUtilities.cc:103 (set (reg/v/f:DI 74 [ perc ])
        (reg:DI 0 ax)) -1 (nil)
    (nil))

(insn 309 308 310 27 ./CppStatUtilities.cc:103 (parallel [
            (set (reg:DI 65 [ prephitmp.1125 ])
                (plus:DI (reg/v/f:DI 74 [ perc ])
                    (const_int 16 [0x10])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(jump_insn 310 309 311 27 ./CppStatUtilities.cc:103 (set (pc)
        (label_ref 215)) -1 (nil)
    (nil))
;; End of basic block 27, registers live:
 (nil)

(barrier 311 310 312)

(note 312 311 313 NOTE_INSN_FUNCTION_END)

(note 313 312 317 ("./CppStatUtilities.cc") 130)

;; Start of basic block 28, registers live: (nil)
(code_label 317 313 320 28 1056 "" [1 uses])

(note 320 317 0 28 [bb 28] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 28, registers live:
 (nil)

