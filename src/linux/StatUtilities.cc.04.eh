
;; Function double scaledL1dist(double*, double*, int) (_Z12scaledL1distPdS_i)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Deleted label in block 5.
Merged 6 and 7 without moving.
Merged 6 and 8 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 14 ("./StatUtilities.cc") 65)

;; Start of basic block 0, registers live: (nil)
(note 14 1 9 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 9 14 10 0 ./StatUtilities.cc:65 (set (reg/v/f:DI 67 [ v1 ])
        (reg:DI 5 di [ v1 ])) -1 (nil)
    (nil))

(insn 10 9 11 0 ./StatUtilities.cc:65 (set (reg/v/f:DI 68 [ v2 ])
        (reg:DI 4 si [ v2 ])) -1 (nil)
    (nil))

(insn 11 10 12 0 ./StatUtilities.cc:65 (set (reg/v:SI 69 [ n ])
        (reg:SI 1 dx [ n ])) -1 (nil)
    (nil))

(note 12 11 16 0 NOTE_INSN_FUNCTION_BEG)

(note 16 12 17 0 ("./StatUtilities.cc") 69)

(insn 17 16 18 0 ./StatUtilities.cc:69 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 69 [ n ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 18 17 20 0 ./StatUtilities.cc:69 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 24)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 20 18 21 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 21 20 22 1 ./StatUtilities.cc:69 (set (reg:DF 63 [ prephitmp.127 ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(jump_insn 22 21 23 1 ./StatUtilities.cc:69 (set (pc)
        (label_ref 50)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 23 22 24)

;; Start of basic block 2, registers live: (nil)
(code_label 24 23 25 2 2 "" [1 uses])

(note 25 24 26 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 26 25 27 2 ./StatUtilities.cc:69 (set (reg/v:DF 65 [ delta ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC0") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(insn 27 26 28 2 ./StatUtilities.cc:69 (set (reg/v:SI 64 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 28 27 29 2 ./StatUtilities.cc:69 (set (reg:DI 62 [ ivtmp.130 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 29 28 30 3 5 "" [1 uses])

(note 30 29 31 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 31 30 32 3 ("./StatUtilities.cc") 65)

(insn 32 31 33 3 ./StatUtilities.cc:65 (set (reg:DI 61 [ D.34803 ])
        (reg:DI 62 [ ivtmp.130 ])) -1 (nil)
    (nil))

(note 33 32 34 3 ("./StatUtilities.cc") 70)

(insn 34 33 35 3 ./StatUtilities.cc:70 (set (reg:DF 71)
        (mem:DF (plus:DI (mult:DI (reg:DI 61 [ D.34803 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 67 [ v1 ])) [2 S8 A64])) -1 (nil)
    (nil))

(insn 35 34 36 3 ./StatUtilities.cc:70 (set (reg:DF 70)
        (minus:DF (reg:DF 71)
            (mem:DF (plus:DI (mult:DI (reg:DI 61 [ D.34803 ])
                        (const_int 8 [0x8]))
                    (reg/v/f:DI 68 [ v2 ])) [2 S8 A64]))) -1 (nil)
    (nil))

(insn 36 35 37 3 ./StatUtilities.cc:70 (set (reg:V2DF 73)
        (mem/u/c/i:V2DF (symbol_ref/u:DI ("*.LC1") [flags 0x2]) [2 S16 A128])) -1 (nil)
    (expr_list:REG_EQUAL (const_vector:V2DF [
                (const_double:DF +NaN [+NaN])
                (const_double:DF 0.0 [0x0.0p+0])
            ])
        (nil)))

(insn 37 36 38 3 ./StatUtilities.cc:70 (parallel [
            (set (reg:DF 72)
                (abs:DF (reg:DF 70)))
            (use (reg:V2DF 73))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (expr_list:REG_EQUAL (abs:DF (reg:DF 70))
        (nil)))

(insn 38 37 39 3 ./StatUtilities.cc:70 (set (reg/v:DF 65 [ delta ])
        (plus:DF (reg/v:DF 65 [ delta ])
            (reg:DF 72))) -1 (nil)
    (nil))

(note 39 38 40 3 ("./StatUtilities.cc") 69)

(insn 40 39 41 3 ./StatUtilities.cc:69 (parallel [
            (set (reg/v:SI 64 [ i ])
                (plus:SI (reg/v:SI 64 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 41 40 42 3 ./StatUtilities.cc:69 (parallel [
            (set (reg:DI 62 [ ivtmp.130 ])
                (plus:DI (reg:DI 62 [ ivtmp.130 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 42 41 43 3 ./StatUtilities.cc:69 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 64 [ i ])
            (reg/v:SI 69 [ n ]))) -1 (nil)
    (nil))

(jump_insn 43 42 45 3 ./StatUtilities.cc:69 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 29)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(note 45 43 46 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 46 45 47 4 ./StatUtilities.cc:69 (set (reg:DF 74)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC2") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 1.0e+3 [0x0.fap+10])
        (nil)))

(insn 47 46 48 4 ./StatUtilities.cc:69 (set (reg:DF 75)
        (mult:DF (reg/v:DF 65 [ delta ])
            (reg:DF 74))) -1 (nil)
    (nil))

(insn 48 47 49 4 ./StatUtilities.cc:69 (set (reg:DF 76)
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC3") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 6.25e-2 [0x0.8p-3])
        (nil)))

(insn 49 48 50 4 ./StatUtilities.cc:69 (set (reg:DF 63 [ prephitmp.127 ])
        (mult:DF (reg:DF 75)
            (reg:DF 76))) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

;; Start of basic block 5, registers live: (nil)
(code_label 50 49 51 5 4 "" [1 uses])

(note 51 50 52 5 [bb 5] NOTE_INSN_BASIC_BLOCK)

(insn 52 51 55 5 ./StatUtilities.cc:69 (set (reg:DF 66 [ <result> ])
        (reg:DF 63 [ prephitmp.127 ])) -1 (nil)
    (nil))

(note 55 52 56 5 NOTE_INSN_FUNCTION_END)

(note 56 55 58 5 ("./StatUtilities.cc") 73)

(insn 58 56 64 5 ./StatUtilities.cc:73 (set (reg/i:DF 21 xmm0)
        (reg:DF 66 [ <result> ])) -1 (nil)
    (nil))

(insn 64 58 0 5 ./StatUtilities.cc:73 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 5, registers live:
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

(note 1 0 7 ("./StatUtilities.cc") 85)

;; Start of basic block 0, registers live: (nil)
(note 7 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 7 4 0 ./StatUtilities.cc:85 (set (reg/v:SI 58 [ __initialize_p ])
        (reg:SI 5 di [ __initialize_p ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./StatUtilities.cc:85 (set (reg/v:SI 59 [ __priority ])
        (reg:SI 4 si [ __priority ])) -1 (nil)
    (nil))

(note 5 4 9 0 NOTE_INSN_FUNCTION_BEG)

(note 9 5 10 0 ("./StatUtilities.cc") 85)

(insn 10 9 11 0 ./StatUtilities.cc:85 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 58 [ __initialize_p ])
            (const_int 1 [0x1]))) -1 (nil)
    (nil))

(jump_insn 11 10 13 0 ./StatUtilities.cc:85 (set (pc)
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

(insn 14 13 15 1 ./StatUtilities.cc:85 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 59 [ __priority ])
            (const_int 65535 [0xffff]))) -1 (nil)
    (nil))

(jump_insn 15 14 17 1 ./StatUtilities.cc:85 (set (pc)
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
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2b096de7aa50 __ioinit>)) -1 (nil)
    (nil))

(call_insn 20 19 21 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitC1Ev") [flags 0x41] <function_decl 0x2b096d635f00 __comp_ctor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 21 20 22 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 1 dx)
        (symbol_ref:DI ("__dso_handle") [flags 0x40] <var_decl 0x2b096e0fff20 __dso_handle>)) -1 (nil)
    (nil))

(insn 22 21 23 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 4 si)
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 23 22 24 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("__tcf_0") [flags 0x3] <function_decl 0x2b096e0e8e00 __tcf_0>)) -1 (nil)
    (nil))

(call_insn/j 24 23 25 2 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:SI 0 ax)
        (call (mem:QI (symbol_ref:DI ("__cxa_atexit") [flags 0x41] <function_decl 0x2b096e0e8f00 __cxa_atexit>) [0 S1 A8])
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

(note 29 28 33 ("./StatUtilities.cc") 85)

;; Start of basic block 3, registers live: (nil)
(code_label 33 29 36 3 47 "" [2 uses])

(note 36 33 0 3 [bb 3] NOTE_INSN_BASIC_BLOCK)
;; End of basic block 3, registers live:
 (nil)


;; Function (static initializers for ./StatUtilities.cc) (_GLOBAL__I__Z3covPdS_i)



try_optimize_cfg iteration 1

Deleting fallthru block 0.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 3 ("./StatUtilities.cc") 86)

(note 3 1 6 0 NOTE_INSN_FUNCTION_BEG)

;; Start of basic block 0, registers live: (nil)
(note 6 3 7 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(note 7 6 8 0 ("./StatUtilities.cc") 86)

(insn 8 7 9 0 ./StatUtilities.cc:86 (set (reg:SI 4 si)
        (const_int 65535 [0xffff])) -1 (nil)
    (nil))

(insn 9 8 10 0 ./StatUtilities.cc:86 (set (reg:SI 5 di)
        (const_int 1 [0x1])) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 ./StatUtilities.cc:86 (call (mem:QI (symbol_ref:DI ("_Z41__static_initialization_and_destruction_0ii") [flags 0x3] <function_decl 0x2b096e0e8d00 __static_initialization_and_destruction_0>) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:SI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("./StatUtilities.cc") 86)


;; Function void __tcf_0(void*) (__tcf_0)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 6 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

;; Start of basic block 0, registers live: (nil)
(note 6 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 6 4 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg/v/f:DI 58 [ D.34721 ])
        (reg:DI 5 di [ D.34721 ])) -1 (nil)
    (nil))

(note 4 3 8 0 NOTE_INSN_FUNCTION_BEG)

(note 8 4 9 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)

(insn 9 8 10 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (set (reg:DI 5 di)
        (symbol_ref:DI ("_ZSt8__ioinit") [flags 0x2] <var_decl 0x2b096de7aa50 __ioinit>)) -1 (nil)
    (nil))

(call_insn/j 10 9 11 0 /usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream:76 (call (mem:QI (symbol_ref:DI ("_ZNSt8ios_base4InitD1Ev") [flags 0x41] <function_decl 0x2b096d63f100 __comp_dtor >) [0 S1 A8])
        (const_int 0 [0x0])) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))
;; End of basic block 0, registers live:
 (nil)

(barrier 11 10 12)

(note 12 11 13 NOTE_INSN_FUNCTION_END)

(note 13 12 0 ("/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../include/c++/4.1.2/iostream") 76)


;; Function double* doubleCol(double**, int, int) (_Z9doubleColPPdii)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 4 and 5 without moving.
Merged 4 and 6 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 8 ("./StatUtilities.cc") 77)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./StatUtilities.cc:77 (set (reg/v/f:DI 64 [ x ])
        (reg:DI 5 di [ x ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./StatUtilities.cc:77 (set (reg/v:SI 65 [ nRows ])
        (reg:SI 4 si [ nRows ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./StatUtilities.cc:77 (set (reg/v:SI 66 [ colIdx ])
        (reg:SI 1 dx [ colIdx ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./StatUtilities.cc") 79)

(insn 11 10 12 0 ./StatUtilities.cc:79 (set (reg:DI 67)
        (sign_extend:DI (reg/v:SI 65 [ nRows ]))) -1 (nil)
    (nil))

(insn 12 11 13 0 ./StatUtilities.cc:79 (parallel [
            (set (reg:DI 68)
                (ashift:DI (reg:DI 67)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 13 12 14 0 ./StatUtilities.cc:79 (set (reg:DI 5 di)
        (reg:DI 68)) -1 (nil)
    (nil))

(call_insn 14 13 15 0 ./StatUtilities.cc:79 (set (reg:DI 0 ax)
        (call (mem:QI (symbol_ref:DI ("malloc") [flags 0x41] <function_decl 0x2b096b997100 malloc>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (expr_list:REG_EH_REGION (const_int 0 [0x0])
        (nil))
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (nil)))

(insn 15 14 16 0 ./StatUtilities.cc:79 (set (reg/f:DI 69)
        (reg:DI 0 ax)) -1 (nil)
    (expr_list:REG_NOALIAS (reg/f:DI 69)
        (nil)))

(insn 16 15 17 0 ./StatUtilities.cc:79 (set (reg:DI 62 [ D.34656 ])
        (reg/f:DI 69)) -1 (nil)
    (nil))

(insn 17 16 18 0 ./StatUtilities.cc:79 (set (reg/v/f:DI 61 [ col ])
        (reg:DI 62 [ D.34656 ])) -1 (nil)
    (nil))

(note 18 17 19 0 ("./StatUtilities.cc") 81)

(insn 19 18 20 0 ./StatUtilities.cc:81 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 65 [ nRows ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 20 19 22 0 ./StatUtilities.cc:81 (set (pc)
        (if_then_else (le (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 38)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 1000 [0x3e8])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 22 20 23 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 23 22 24 1 ./StatUtilities.cc:81 (set (reg:DI 70)
        (sign_extend:DI (reg/v:SI 66 [ colIdx ]))) -1 (nil)
    (nil))

(insn 24 23 25 1 ./StatUtilities.cc:81 (parallel [
            (set (reg:DI 59 [ pretmp.244 ])
                (ashift:DI (reg:DI 70)
                    (const_int 3 [0x3])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 25 24 26 1 ./StatUtilities.cc:81 (set (reg/v:SI 60 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 26 25 27 1 ./StatUtilities.cc:81 (set (reg:DI 58 [ ivtmp.247 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

;; Start of basic block 2, registers live: (nil)
(code_label 27 26 28 2 56 "" [1 uses])

(note 28 27 29 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(note 29 28 30 2 ("./StatUtilities.cc") 82)

(insn 30 29 31 2 ./StatUtilities.cc:82 (set (reg:DI 71)
        (mem/f:DI (plus:DI (mult:DI (reg:DI 58 [ ivtmp.247 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 64 [ x ])) [5 S8 A64])) -1 (nil)
    (nil))

(insn 31 30 32 2 ./StatUtilities.cc:82 (set (reg:DF 72)
        (mem:DF (plus:DI (reg:DI 71)
                (reg:DI 59 [ pretmp.244 ])) [2 S8 A64])) -1 (nil)
    (nil))

(insn 32 31 33 2 ./StatUtilities.cc:82 (set (mem:DF (plus:DI (mult:DI (reg:DI 58 [ ivtmp.247 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 61 [ col ])) [2 S8 A64])
        (reg:DF 72)) -1 (nil)
    (nil))

(note 33 32 34 2 ("./StatUtilities.cc") 81)

(insn 34 33 35 2 ./StatUtilities.cc:81 (parallel [
            (set (reg/v:SI 60 [ i ])
                (plus:SI (reg/v:SI 60 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 35 34 36 2 ./StatUtilities.cc:81 (parallel [
            (set (reg:DI 58 [ ivtmp.247 ])
                (plus:DI (reg:DI 58 [ ivtmp.247 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 36 35 37 2 ./StatUtilities.cc:81 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 60 [ i ])
            (reg/v:SI 65 [ nRows ]))) -1 (nil)
    (nil))

(jump_insn 37 36 38 2 ./StatUtilities.cc:81 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 27)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 38 37 39 3 54 "" [1 uses])

(note 39 38 40 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(insn 40 39 43 3 ./StatUtilities.cc:81 (set (reg:DI 63 [ <result> ])
        (reg/v/f:DI 61 [ col ])) -1 (nil)
    (nil))

(note 43 40 44 3 NOTE_INSN_FUNCTION_END)

(note 44 43 46 3 ("./StatUtilities.cc") 85)

(insn 46 44 52 3 ./StatUtilities.cc:85 (set (reg/i:DI 0 ax)
        (reg:DI 63 [ <result> ])) -1 (nil)
    (nil))

(insn 52 46 0 3 ./StatUtilities.cc:85 (use (reg/i:DI 0 ax)) -1 (nil)
    (nil))
;; End of basic block 3, registers live:
 (nil)


;; Function double rcov(double*, double*, int) (_Z4rcovPdS_i)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 5 and 6 without moving.
Merged 5 and 7 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 8 ("./StatUtilities.cc") 49)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./StatUtilities.cc:49 (set (reg/v/f:DI 65 [ x ])
        (reg:DI 5 di [ x ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./StatUtilities.cc:49 (set (reg/v/f:DI 66 [ y ])
        (reg:DI 4 si [ y ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./StatUtilities.cc:49 (set (reg/v:SI 67 [ n ])
        (reg:SI 1 dx [ n ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./StatUtilities.cc") 53)

(insn 11 10 12 0 ./StatUtilities.cc:53 (set (reg:SI 4 si)
        (reg/v:SI 67 [ n ])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./StatUtilities.cc:53 (set (reg:DI 5 di)
        (reg/v/f:DI 65 [ x ])) -1 (nil)
    (nil))

(call_insn 13 12 14 0 ./StatUtilities.cc:53 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("_Z6medianPdi") [flags 0x41] <function_decl 0x2b096e01c200 median>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 14 13 15 0 ./StatUtilities.cc:53 (set (reg/v:DF 62 [ mx ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 15 14 16 0 ("./StatUtilities.cc") 54)

(insn 16 15 17 0 ./StatUtilities.cc:54 (set (reg:SI 4 si)
        (reg/v:SI 67 [ n ])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./StatUtilities.cc:54 (set (reg:DI 5 di)
        (reg/v/f:DI 66 [ y ])) -1 (nil)
    (nil))

(call_insn 18 17 19 0 ./StatUtilities.cc:54 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("_Z6medianPdi") [flags 0x41] <function_decl 0x2b096e01c200 median>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 19 18 20 0 ./StatUtilities.cc:54 (set (reg/v:DF 61 [ my ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 20 19 21 0 ("./StatUtilities.cc") 56)

(insn 21 20 22 0 ./StatUtilities.cc:56 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 67 [ n ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 22 21 24 0 ./StatUtilities.cc:56 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 28)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 24 22 25 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 25 24 26 1 ./StatUtilities.cc:56 (set (reg/v:DF 63 [ z ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC4") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(jump_insn 26 25 27 1 ./StatUtilities.cc:56 (set (pc)
        (label_ref 49)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 27 26 28)

;; Start of basic block 2, registers live: (nil)
(code_label 28 27 29 2 92 "" [1 uses])

(note 29 28 30 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 30 29 31 2 ./StatUtilities.cc:56 (set (reg/v:DF 63 [ z ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC4") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(insn 31 30 32 2 ./StatUtilities.cc:56 (set (reg/v:SI 60 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 32 31 33 2 ./StatUtilities.cc:56 (set (reg:DI 59 [ ivtmp.283 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 33 32 34 3 95 "" [1 uses])

(note 34 33 35 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 35 34 36 3 ("./StatUtilities.cc") 49)

(insn 36 35 37 3 ./StatUtilities.cc:49 (set (reg:DI 58 [ D.35050 ])
        (reg:DI 59 [ ivtmp.283 ])) -1 (nil)
    (nil))

(note 37 36 38 3 ("./StatUtilities.cc") 57)

(insn 38 37 39 3 ./StatUtilities.cc:57 (set (reg:DF 69)
        (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.35050 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 65 [ x ])) [2 S8 A64])) -1 (nil)
    (nil))

(insn 39 38 40 3 ./StatUtilities.cc:57 (set (reg:DF 68)
        (minus:DF (reg:DF 69)
            (reg/v:DF 62 [ mx ]))) -1 (nil)
    (nil))

(insn 40 39 41 3 ./StatUtilities.cc:57 (set (reg:DF 71)
        (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.35050 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 66 [ y ])) [2 S8 A64])) -1 (nil)
    (nil))

(insn 41 40 42 3 ./StatUtilities.cc:57 (set (reg:DF 70)
        (minus:DF (reg:DF 71)
            (reg/v:DF 61 [ my ]))) -1 (nil)
    (nil))

(insn 42 41 43 3 ./StatUtilities.cc:57 (set (reg:DF 72)
        (mult:DF (reg:DF 68)
            (reg:DF 70))) -1 (nil)
    (nil))

(insn 43 42 44 3 ./StatUtilities.cc:57 (set (reg/v:DF 63 [ z ])
        (plus:DF (reg/v:DF 63 [ z ])
            (reg:DF 72))) -1 (nil)
    (nil))

(note 44 43 45 3 ("./StatUtilities.cc") 56)

(insn 45 44 46 3 ./StatUtilities.cc:56 (parallel [
            (set (reg/v:SI 60 [ i ])
                (plus:SI (reg/v:SI 60 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 46 45 47 3 ./StatUtilities.cc:56 (parallel [
            (set (reg:DI 59 [ ivtmp.283 ])
                (plus:DI (reg:DI 59 [ ivtmp.283 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 47 46 48 3 ./StatUtilities.cc:56 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 60 [ i ])
            (reg/v:SI 67 [ n ]))) -1 (nil)
    (nil))

(jump_insn 48 47 49 3 ./StatUtilities.cc:56 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 33)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 49 48 50 4 94 "" [1 uses])

(note 50 49 51 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 51 50 52 4 ./StatUtilities.cc:56 (parallel [
            (set (reg:SI 74)
                (plus:SI (reg/v:SI 67 [ n ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 52 51 53 4 ./StatUtilities.cc:56 (set (reg:DF 75)
        (float:DF (reg:SI 74))) -1 (nil)
    (nil))

(insn 53 52 54 4 ./StatUtilities.cc:56 (set (reg:DF 73)
        (div:DF (reg/v:DF 63 [ z ])
            (reg:DF 75))) -1 (nil)
    (nil))

(insn 54 53 57 4 ./StatUtilities.cc:56 (set (reg:DF 64 [ <result> ])
        (reg:DF 73)) -1 (nil)
    (nil))

(note 57 54 58 4 NOTE_INSN_FUNCTION_END)

(note 58 57 60 4 ("./StatUtilities.cc") 60)

(insn 60 58 66 4 ./StatUtilities.cc:60 (set (reg/i:DF 21 xmm0)
        (reg:DF 64 [ <result> ])) -1 (nil)
    (nil))

(insn 66 60 0 4 ./StatUtilities.cc:60 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)


;; Function double cov(double*, double*, int) (_Z3covPdS_i)



try_optimize_cfg iteration 1

Merged 0 and 1 without moving.
Deleted label in block 2.
Merged 5 and 6 without moving.
Merged 5 and 7 without moving.


try_optimize_cfg iteration 2



try_optimize_cfg iteration 1

(note 1 0 8 ("./StatUtilities.cc") 34)

;; Start of basic block 0, registers live: (nil)
(note 8 1 3 0 [bb 0] NOTE_INSN_BASIC_BLOCK)

(insn 3 8 4 0 ./StatUtilities.cc:34 (set (reg/v/f:DI 65 [ x ])
        (reg:DI 5 di [ x ])) -1 (nil)
    (nil))

(insn 4 3 5 0 ./StatUtilities.cc:34 (set (reg/v/f:DI 66 [ y ])
        (reg:DI 4 si [ y ])) -1 (nil)
    (nil))

(insn 5 4 6 0 ./StatUtilities.cc:34 (set (reg/v:SI 67 [ n ])
        (reg:SI 1 dx [ n ])) -1 (nil)
    (nil))

(note 6 5 10 0 NOTE_INSN_FUNCTION_BEG)

(note 10 6 11 0 ("./StatUtilities.cc") 38)

(insn 11 10 12 0 ./StatUtilities.cc:38 (set (reg:SI 4 si)
        (reg/v:SI 67 [ n ])) -1 (nil)
    (nil))

(insn 12 11 13 0 ./StatUtilities.cc:38 (set (reg:DI 5 di)
        (reg/v/f:DI 65 [ x ])) -1 (nil)
    (nil))

(call_insn 13 12 14 0 ./StatUtilities.cc:38 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("mean") [flags 0x41] <function_decl 0x2b096de7d700 mean>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 14 13 15 0 ./StatUtilities.cc:38 (set (reg/v:DF 62 [ mx ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 15 14 16 0 ("./StatUtilities.cc") 39)

(insn 16 15 17 0 ./StatUtilities.cc:39 (set (reg:SI 4 si)
        (reg/v:SI 67 [ n ])) -1 (nil)
    (nil))

(insn 17 16 18 0 ./StatUtilities.cc:39 (set (reg:DI 5 di)
        (reg/v/f:DI 66 [ y ])) -1 (nil)
    (nil))

(call_insn 18 17 19 0 ./StatUtilities.cc:39 (set (reg:DF 21 xmm0)
        (call (mem:QI (symbol_ref:DI ("mean") [flags 0x41] <function_decl 0x2b096de7d700 mean>) [0 S1 A8])
            (const_int 0 [0x0]))) -1 (nil)
    (nil)
    (expr_list:REG_DEP_TRUE (use (reg:DI 5 di))
        (expr_list:REG_DEP_TRUE (use (reg:SI 4 si))
            (nil))))

(insn 19 18 20 0 ./StatUtilities.cc:39 (set (reg/v:DF 61 [ my ])
        (reg:DF 21 xmm0)) -1 (nil)
    (nil))

(note 20 19 21 0 ("./StatUtilities.cc") 41)

(insn 21 20 22 0 ./StatUtilities.cc:41 (set (reg:CCNO 17 flags)
        (compare:CCNO (reg/v:SI 67 [ n ])
            (const_int 0 [0x0]))) -1 (nil)
    (nil))

(jump_insn 22 21 24 0 ./StatUtilities.cc:41 (set (pc)
        (if_then_else (gt (reg:CCNO 17 flags)
                (const_int 0 [0x0]))
            (label_ref 28)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 0, registers live:
 (nil)

;; Start of basic block 1, registers live: (nil)
(note 24 22 25 1 [bb 1] NOTE_INSN_BASIC_BLOCK)

(insn 25 24 26 1 ./StatUtilities.cc:41 (set (reg/v:DF 63 [ z ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC5") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(jump_insn 26 25 27 1 ./StatUtilities.cc:41 (set (pc)
        (label_ref 49)) -1 (nil)
    (nil))
;; End of basic block 1, registers live:
 (nil)

(barrier 27 26 28)

;; Start of basic block 2, registers live: (nil)
(code_label 28 27 29 2 131 "" [1 uses])

(note 29 28 30 2 [bb 2] NOTE_INSN_BASIC_BLOCK)

(insn 30 29 31 2 ./StatUtilities.cc:41 (set (reg/v:DF 63 [ z ])
        (mem/u/c/i:DF (symbol_ref/u:DI ("*.LC5") [flags 0x2]) [2 S8 A64])) -1 (nil)
    (expr_list:REG_EQUAL (const_double:DF 0.0 [0x0.0p+0])
        (nil)))

(insn 31 30 32 2 ./StatUtilities.cc:41 (set (reg/v:SI 60 [ i ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))

(insn 32 31 33 2 ./StatUtilities.cc:41 (set (reg:DI 59 [ ivtmp.319 ])
        (const_int 0 [0x0])) -1 (nil)
    (nil))
;; End of basic block 2, registers live:
 (nil)

;; Start of basic block 3, registers live: (nil)
(code_label 33 32 34 3 134 "" [1 uses])

(note 34 33 35 3 [bb 3] NOTE_INSN_BASIC_BLOCK)

(note 35 34 36 3 ("./StatUtilities.cc") 34)

(insn 36 35 37 3 ./StatUtilities.cc:34 (set (reg:DI 58 [ D.35127 ])
        (reg:DI 59 [ ivtmp.319 ])) -1 (nil)
    (nil))

(note 37 36 38 3 ("./StatUtilities.cc") 42)

(insn 38 37 39 3 ./StatUtilities.cc:42 (set (reg:DF 69)
        (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.35127 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 65 [ x ])) [2 S8 A64])) -1 (nil)
    (nil))

(insn 39 38 40 3 ./StatUtilities.cc:42 (set (reg:DF 68)
        (minus:DF (reg:DF 69)
            (reg/v:DF 62 [ mx ]))) -1 (nil)
    (nil))

(insn 40 39 41 3 ./StatUtilities.cc:42 (set (reg:DF 71)
        (mem:DF (plus:DI (mult:DI (reg:DI 58 [ D.35127 ])
                    (const_int 8 [0x8]))
                (reg/v/f:DI 66 [ y ])) [2 S8 A64])) -1 (nil)
    (nil))

(insn 41 40 42 3 ./StatUtilities.cc:42 (set (reg:DF 70)
        (minus:DF (reg:DF 71)
            (reg/v:DF 61 [ my ]))) -1 (nil)
    (nil))

(insn 42 41 43 3 ./StatUtilities.cc:42 (set (reg:DF 72)
        (mult:DF (reg:DF 68)
            (reg:DF 70))) -1 (nil)
    (nil))

(insn 43 42 44 3 ./StatUtilities.cc:42 (set (reg/v:DF 63 [ z ])
        (plus:DF (reg/v:DF 63 [ z ])
            (reg:DF 72))) -1 (nil)
    (nil))

(note 44 43 45 3 ("./StatUtilities.cc") 41)

(insn 45 44 46 3 ./StatUtilities.cc:41 (parallel [
            (set (reg/v:SI 60 [ i ])
                (plus:SI (reg/v:SI 60 [ i ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 46 45 47 3 ./StatUtilities.cc:41 (parallel [
            (set (reg:DI 59 [ ivtmp.319 ])
                (plus:DI (reg:DI 59 [ ivtmp.319 ])
                    (const_int 1 [0x1])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 47 46 48 3 ./StatUtilities.cc:41 (set (reg:CCZ 17 flags)
        (compare:CCZ (reg/v:SI 60 [ i ])
            (reg/v:SI 67 [ n ]))) -1 (nil)
    (nil))

(jump_insn 48 47 49 3 ./StatUtilities.cc:41 (set (pc)
        (if_then_else (ne (reg:CCZ 17 flags)
                (const_int 0 [0x0]))
            (label_ref 33)
            (pc))) -1 (nil)
    (expr_list:REG_BR_PROB (const_int 9000 [0x2328])
        (nil)))
;; End of basic block 3, registers live:
 (nil)

;; Start of basic block 4, registers live: (nil)
(code_label 49 48 50 4 133 "" [1 uses])

(note 50 49 51 4 [bb 4] NOTE_INSN_BASIC_BLOCK)

(insn 51 50 52 4 ./StatUtilities.cc:41 (parallel [
            (set (reg:SI 74)
                (plus:SI (reg/v:SI 67 [ n ])
                    (const_int -1 [0xffffffffffffffff])))
            (clobber (reg:CC 17 flags))
        ]) -1 (nil)
    (nil))

(insn 52 51 53 4 ./StatUtilities.cc:41 (set (reg:DF 75)
        (float:DF (reg:SI 74))) -1 (nil)
    (nil))

(insn 53 52 54 4 ./StatUtilities.cc:41 (set (reg:DF 73)
        (div:DF (reg/v:DF 63 [ z ])
            (reg:DF 75))) -1 (nil)
    (nil))

(insn 54 53 57 4 ./StatUtilities.cc:41 (set (reg:DF 64 [ <result> ])
        (reg:DF 73)) -1 (nil)
    (nil))

(note 57 54 58 4 NOTE_INSN_FUNCTION_END)

(note 58 57 60 4 ("./StatUtilities.cc") 45)

(insn 60 58 66 4 ./StatUtilities.cc:45 (set (reg/i:DF 21 xmm0)
        (reg:DF 64 [ <result> ])) -1 (nil)
    (nil))

(insn 66 60 0 4 ./StatUtilities.cc:45 (use (reg/i:DF 21 xmm0)) -1 (nil)
    (nil))
;; End of basic block 4, registers live:
 (nil)

