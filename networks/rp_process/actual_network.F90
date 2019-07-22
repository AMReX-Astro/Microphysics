module actual_network

  use physical_constants, only: ERG_PER_MeV
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

  real(rt), parameter :: avo = 6.0221417930d23
  real(rt), parameter :: c_light = 2.99792458d10
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740d-12
  real(rt), parameter :: mev2erg = ev2erg*1.0d6
  real(rt), parameter :: mev2gr  = mev2erg/c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184d-24
  real(rt), parameter :: mass_proton   = 1.67262163783d-24
  real(rt), parameter :: mass_electron = 9.10938215450d-28

  integer, parameter :: nrates = 1174
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 197
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 197

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 1174
  integer, parameter :: number_reaclib_sets = 1601

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jd   = 2
  integer, parameter :: jhe3   = 3
  integer, parameter :: jhe4   = 4
  integer, parameter :: jli7   = 5
  integer, parameter :: jbe7   = 6
  integer, parameter :: jbe8   = 7
  integer, parameter :: jb8   = 8
  integer, parameter :: jc12   = 9
  integer, parameter :: jc13   = 10
  integer, parameter :: jn12   = 11
  integer, parameter :: jn13   = 12
  integer, parameter :: jn14   = 13
  integer, parameter :: jn15   = 14
  integer, parameter :: jo14   = 15
  integer, parameter :: jo15   = 16
  integer, parameter :: jo16   = 17
  integer, parameter :: jo17   = 18
  integer, parameter :: jo18   = 19
  integer, parameter :: jf16   = 20
  integer, parameter :: jf17   = 21
  integer, parameter :: jf18   = 22
  integer, parameter :: jf19   = 23
  integer, parameter :: jf20   = 24
  integer, parameter :: jne17   = 25
  integer, parameter :: jne18   = 26
  integer, parameter :: jne19   = 27
  integer, parameter :: jne20   = 28
  integer, parameter :: jne21   = 29
  integer, parameter :: jne22   = 30
  integer, parameter :: jna19   = 31
  integer, parameter :: jna20   = 32
  integer, parameter :: jna21   = 33
  integer, parameter :: jna22   = 34
  integer, parameter :: jna23   = 35
  integer, parameter :: jna24   = 36
  integer, parameter :: jmg21   = 37
  integer, parameter :: jmg22   = 38
  integer, parameter :: jmg23   = 39
  integer, parameter :: jmg24   = 40
  integer, parameter :: jmg25   = 41
  integer, parameter :: jmg26   = 42
  integer, parameter :: jal22   = 43
  integer, parameter :: jal23   = 44
  integer, parameter :: jal24   = 45
  integer, parameter :: jal25   = 46
  integer, parameter :: jal26   = 47
  integer, parameter :: jal27   = 48
  integer, parameter :: jal28   = 49
  integer, parameter :: jsi24   = 50
  integer, parameter :: jsi25   = 51
  integer, parameter :: jsi26   = 52
  integer, parameter :: jsi27   = 53
  integer, parameter :: jsi28   = 54
  integer, parameter :: jsi29   = 55
  integer, parameter :: jsi30   = 56
  integer, parameter :: jp26   = 57
  integer, parameter :: jp27   = 58
  integer, parameter :: jp28   = 59
  integer, parameter :: jp29   = 60
  integer, parameter :: jp30   = 61
  integer, parameter :: jp31   = 62
  integer, parameter :: jp32   = 63
  integer, parameter :: jp33   = 64
  integer, parameter :: js28   = 65
  integer, parameter :: js29   = 66
  integer, parameter :: js30   = 67
  integer, parameter :: js31   = 68
  integer, parameter :: js32   = 69
  integer, parameter :: js33   = 70
  integer, parameter :: js34   = 71
  integer, parameter :: js35   = 72
  integer, parameter :: jcl29   = 73
  integer, parameter :: jcl30   = 74
  integer, parameter :: jcl31   = 75
  integer, parameter :: jcl32   = 76
  integer, parameter :: jcl33   = 77
  integer, parameter :: jcl34   = 78
  integer, parameter :: jcl35   = 79
  integer, parameter :: jcl36   = 80
  integer, parameter :: jcl37   = 81
  integer, parameter :: jar31   = 82
  integer, parameter :: jar32   = 83
  integer, parameter :: jar33   = 84
  integer, parameter :: jar34   = 85
  integer, parameter :: jar35   = 86
  integer, parameter :: jar36   = 87
  integer, parameter :: jar37   = 88
  integer, parameter :: jar38   = 89
  integer, parameter :: jar39   = 90
  integer, parameter :: jk33   = 91
  integer, parameter :: jk34   = 92
  integer, parameter :: jk35   = 93
  integer, parameter :: jk36   = 94
  integer, parameter :: jk37   = 95
  integer, parameter :: jk38   = 96
  integer, parameter :: jk39   = 97
  integer, parameter :: jk40   = 98
  integer, parameter :: jk41   = 99
  integer, parameter :: jca34   = 100
  integer, parameter :: jca35   = 101
  integer, parameter :: jca36   = 102
  integer, parameter :: jca37   = 103
  integer, parameter :: jca38   = 104
  integer, parameter :: jca39   = 105
  integer, parameter :: jca40   = 106
  integer, parameter :: jca41   = 107
  integer, parameter :: jca42   = 108
  integer, parameter :: jca43   = 109
  integer, parameter :: jca44   = 110
  integer, parameter :: jsc36   = 111
  integer, parameter :: jsc37   = 112
  integer, parameter :: jsc38   = 113
  integer, parameter :: jsc39   = 114
  integer, parameter :: jsc40   = 115
  integer, parameter :: jsc41   = 116
  integer, parameter :: jsc42   = 117
  integer, parameter :: jsc43   = 118
  integer, parameter :: jsc44   = 119
  integer, parameter :: jsc45   = 120
  integer, parameter :: jsc46   = 121
  integer, parameter :: jti38   = 122
  integer, parameter :: jti39   = 123
  integer, parameter :: jti40   = 124
  integer, parameter :: jti41   = 125
  integer, parameter :: jti42   = 126
  integer, parameter :: jti43   = 127
  integer, parameter :: jti44   = 128
  integer, parameter :: jti45   = 129
  integer, parameter :: jti46   = 130
  integer, parameter :: jti47   = 131
  integer, parameter :: jti48   = 132
  integer, parameter :: jv40   = 133
  integer, parameter :: jv41   = 134
  integer, parameter :: jv42   = 135
  integer, parameter :: jv43   = 136
  integer, parameter :: jv44   = 137
  integer, parameter :: jv45   = 138
  integer, parameter :: jv46   = 139
  integer, parameter :: jv47   = 140
  integer, parameter :: jv48   = 141
  integer, parameter :: jv49   = 142
  integer, parameter :: jv50   = 143
  integer, parameter :: jcr42   = 144
  integer, parameter :: jcr43   = 145
  integer, parameter :: jcr44   = 146
  integer, parameter :: jcr45   = 147
  integer, parameter :: jcr46   = 148
  integer, parameter :: jcr47   = 149
  integer, parameter :: jcr48   = 150
  integer, parameter :: jcr49   = 151
  integer, parameter :: jcr50   = 152
  integer, parameter :: jcr51   = 153
  integer, parameter :: jcr52   = 154
  integer, parameter :: jmn44   = 155
  integer, parameter :: jmn45   = 156
  integer, parameter :: jmn46   = 157
  integer, parameter :: jmn47   = 158
  integer, parameter :: jmn48   = 159
  integer, parameter :: jmn49   = 160
  integer, parameter :: jmn50   = 161
  integer, parameter :: jmn51   = 162
  integer, parameter :: jmn52   = 163
  integer, parameter :: jmn53   = 164
  integer, parameter :: jmn54   = 165
  integer, parameter :: jmn55   = 166
  integer, parameter :: jfe45   = 167
  integer, parameter :: jfe46   = 168
  integer, parameter :: jfe47   = 169
  integer, parameter :: jfe48   = 170
  integer, parameter :: jfe49   = 171
  integer, parameter :: jfe50   = 172
  integer, parameter :: jfe51   = 173
  integer, parameter :: jfe52   = 174
  integer, parameter :: jfe53   = 175
  integer, parameter :: jfe54   = 176
  integer, parameter :: jfe55   = 177
  integer, parameter :: jfe56   = 178
  integer, parameter :: jco47   = 179
  integer, parameter :: jco48   = 180
  integer, parameter :: jco49   = 181
  integer, parameter :: jco50   = 182
  integer, parameter :: jco51   = 183
  integer, parameter :: jco52   = 184
  integer, parameter :: jco53   = 185
  integer, parameter :: jco54   = 186
  integer, parameter :: jco55   = 187
  integer, parameter :: jco56   = 188
  integer, parameter :: jni48   = 189
  integer, parameter :: jni49   = 190
  integer, parameter :: jni50   = 191
  integer, parameter :: jni51   = 192
  integer, parameter :: jni52   = 193
  integer, parameter :: jni53   = 194
  integer, parameter :: jni54   = 195
  integer, parameter :: jni55   = 196
  integer, parameter :: jni56   = 197

  ! Reactions
  integer, parameter :: k_be7__li7__weak__electron_capture   = 1
  integer, parameter :: k_b8__be8__weak__wc17   = 2
  integer, parameter :: k_o14__n14__weak__wc12   = 3
  integer, parameter :: k_o15__n15__weak__wc12   = 4
  integer, parameter :: k_f17__o17__weak__wc12   = 5
  integer, parameter :: k_f18__o18__weak__wc12   = 6
  integer, parameter :: k_f20__ne20__weak__wc12   = 7
  integer, parameter :: k_ne18__f18__weak__wc12   = 8
  integer, parameter :: k_ne19__f19__weak__wc12   = 9
  integer, parameter :: k_he3__p_d   = 10
  integer, parameter :: k_he4__d_d   = 11
  integer, parameter :: k_be7__he4_he3   = 12
  integer, parameter :: k_b8__p_be7   = 13
  integer, parameter :: k_b8__he4_he4__weak__wc12   = 14
  integer, parameter :: k_n13__p_c12   = 15
  integer, parameter :: k_o14__p_n13   = 16
  integer, parameter :: k_o15__p_n14   = 17
  integer, parameter :: k_o16__p_n15   = 18
  integer, parameter :: k_o16__he4_c12   = 19
  integer, parameter :: k_f17__p_o16   = 20
  integer, parameter :: k_f18__p_o17   = 21
  integer, parameter :: k_f18__he4_n14   = 22
  integer, parameter :: k_f19__p_o18   = 23
  integer, parameter :: k_f19__he4_n15   = 24
  integer, parameter :: k_ne18__p_f17   = 25
  integer, parameter :: k_ne18__he4_o14   = 26
  integer, parameter :: k_ne19__p_f18   = 27
  integer, parameter :: k_ne19__he4_o15   = 28
  integer, parameter :: k_ne20__p_f19   = 29
  integer, parameter :: k_ne20__he4_o16   = 30
  integer, parameter :: k_ne21__p_f20   = 31
  integer, parameter :: k_ne21__he4_o17   = 32
  integer, parameter :: k_c12__he4_he4_he4   = 33
  integer, parameter :: k_p_p__d__weak__bet_pos_   = 34
  integer, parameter :: k_p_p__d__weak__electron_capture   = 35
  integer, parameter :: k_p_d__he3   = 36
  integer, parameter :: k_d_d__he4   = 37
  integer, parameter :: k_p_he3__he4__weak__bet_pos_   = 38
  integer, parameter :: k_he4_he3__be7   = 39
  integer, parameter :: k_p_be7__b8   = 40
  integer, parameter :: k_p_c12__n13   = 41
  integer, parameter :: k_he4_c12__o16   = 42
  integer, parameter :: k_p_n13__o14   = 43
  integer, parameter :: k_p_n14__o15   = 44
  integer, parameter :: k_he4_n14__f18   = 45
  integer, parameter :: k_p_n15__o16   = 46
  integer, parameter :: k_he4_n15__f19   = 47
  integer, parameter :: k_he4_o14__ne18   = 48
  integer, parameter :: k_he4_o15__ne19   = 49
  integer, parameter :: k_p_o16__f17   = 50
  integer, parameter :: k_he4_o16__ne20   = 51
  integer, parameter :: k_p_o17__f18   = 52
  integer, parameter :: k_he4_o17__ne21   = 53
  integer, parameter :: k_p_o18__f19   = 54
  integer, parameter :: k_p_f17__ne18   = 55
  integer, parameter :: k_p_f18__ne19   = 56
  integer, parameter :: k_p_f19__ne20   = 57
  integer, parameter :: k_p_f20__ne21   = 58
  integer, parameter :: k_d_he3__p_he4   = 59
  integer, parameter :: k_p_he4__d_he3   = 60
  integer, parameter :: k_he4_he4__p_li7   = 61
  integer, parameter :: k_p_li7__he4_he4   = 62
  integer, parameter :: k_he4_c12__p_n15   = 63
  integer, parameter :: k_c12_c12__he4_ne20   = 64
  integer, parameter :: k_he4_n13__p_o16   = 65
  integer, parameter :: k_he4_n14__p_o17   = 66
  integer, parameter :: k_p_n15__he4_c12   = 67
  integer, parameter :: k_he4_n15__p_o18   = 68
  integer, parameter :: k_he4_o14__p_f17   = 69
  integer, parameter :: k_he4_o15__p_f18   = 70
  integer, parameter :: k_p_o16__he4_n13   = 71
  integer, parameter :: k_he4_o16__p_f19   = 72
  integer, parameter :: k_p_o17__he4_n14   = 73
  integer, parameter :: k_he4_o17__p_f20   = 74
  integer, parameter :: k_p_o18__he4_n15   = 75
  integer, parameter :: k_p_f17__he4_o14   = 76
  integer, parameter :: k_he4_f17__p_ne20   = 77
  integer, parameter :: k_p_f18__he4_o15   = 78
  integer, parameter :: k_he4_f18__p_ne21   = 79
  integer, parameter :: k_p_f19__he4_o16   = 80
  integer, parameter :: k_p_f20__he4_o17   = 81
  integer, parameter :: k_p_ne20__he4_f17   = 82
  integer, parameter :: k_he4_ne20__c12_c12   = 83
  integer, parameter :: k_p_ne21__he4_f18   = 84
  integer, parameter :: k_he3_he3__p_p_he4   = 85
  integer, parameter :: k_d_be7__p_he4_he4   = 86
  integer, parameter :: k_he3_be7__p_p_he4_he4   = 87
  integer, parameter :: k_he4_he4_he4__c12   = 88
  integer, parameter :: k_p_p_he4__he3_he3   = 89
  integer, parameter :: k_p_he4_he4__d_be7   = 90
  integer, parameter :: k_p_p_he4_he4__he3_be7   = 91
  integer, parameter :: k_n13__c13__weak__wc12   = 92
  integer, parameter :: k_n14__p_c13   = 93
  integer, parameter :: k_p_o15__f16   = 94
  integer, parameter :: k_p_o15__he4_n12   = 95
  integer, parameter :: k_he4_o18__ne22   = 96
  integer, parameter :: k_he4_f17__na21   = 97
  integer, parameter :: k_he4_f18__na22   = 98
  integer, parameter :: k_he4_f19__na23   = 99
  integer, parameter :: k_he4_f19__p_ne22   = 100
  integer, parameter :: k_he4_f20__na24   = 101
  integer, parameter :: k_p_ne18__na19   = 102
  integer, parameter :: k_he4_ne18__mg22   = 103
  integer, parameter :: k_he4_ne18__p_na21   = 104
  integer, parameter :: k_p_ne19__na20   = 105
  integer, parameter :: k_he4_ne19__mg23   = 106
  integer, parameter :: k_p_ne19__he4_f16   = 107
  integer, parameter :: k_he4_ne19__p_na22   = 108
  integer, parameter :: k_p_ne20__na21   = 109
  integer, parameter :: k_he4_ne20__mg24   = 110
  integer, parameter :: k_he4_ne20__p_na23   = 111
  integer, parameter :: k_p_ne21__na22   = 112
  integer, parameter :: k_he4_ne21__mg25   = 113
  integer, parameter :: k_he4_ne21__p_na24   = 114
  integer, parameter :: k_p_c13__n14   = 115
  integer, parameter :: k_n12__c12__weak__wc12   = 116
  integer, parameter :: k_he4_n12__p_o15   = 117
  integer, parameter :: k_f16__p_o15   = 118
  integer, parameter :: k_he4_f16__na20   = 119
  integer, parameter :: k_he4_f16__p_ne19   = 120
  integer, parameter :: k_p_ne22__na23   = 121
  integer, parameter :: k_he4_ne22__mg26   = 122
  integer, parameter :: k_p_ne22__he4_f19   = 123
  integer, parameter :: k_na21__ne21__weak__wc12   = 124
  integer, parameter :: k_na21__p_ne20   = 125
  integer, parameter :: k_na21__he4_f17   = 126
  integer, parameter :: k_p_na21__mg22   = 127
  integer, parameter :: k_he4_na21__al25   = 128
  integer, parameter :: k_p_na21__he4_ne18   = 129
  integer, parameter :: k_he4_na21__p_mg24   = 130
  integer, parameter :: k_na22__ne22__weak__wc12   = 131
  integer, parameter :: k_na22__p_ne21   = 132
  integer, parameter :: k_na22__he4_f18   = 133
  integer, parameter :: k_p_na22__mg23   = 134
  integer, parameter :: k_he4_na22__al26   = 135
  integer, parameter :: k_p_na22__he4_ne19   = 136
  integer, parameter :: k_he4_na22__p_mg25   = 137
  integer, parameter :: k_na23__p_ne22   = 138
  integer, parameter :: k_na23__he4_f19   = 139
  integer, parameter :: k_p_na23__mg24   = 140
  integer, parameter :: k_he4_na23__al27   = 141
  integer, parameter :: k_p_na23__he4_ne20   = 142
  integer, parameter :: k_he4_na23__p_mg26   = 143
  integer, parameter :: k_p_na24__mg25   = 144
  integer, parameter :: k_he4_na24__al28   = 145
  integer, parameter :: k_p_na24__he4_ne21   = 146
  integer, parameter :: k_na19__ne19__weak__bqa_pos_   = 147
  integer, parameter :: k_na19__p_ne18   = 148
  integer, parameter :: k_he4_na19__al23   = 149
  integer, parameter :: k_he4_na19__p_mg22   = 150
  integer, parameter :: k_mg22__na22__weak__wc12   = 151
  integer, parameter :: k_mg22__p_na21   = 152
  integer, parameter :: k_mg22__he4_ne18   = 153
  integer, parameter :: k_p_mg22__al23   = 154
  integer, parameter :: k_he4_mg22__si26   = 155
  integer, parameter :: k_p_mg22__he4_na19   = 156
  integer, parameter :: k_he4_mg22__p_al25   = 157
  integer, parameter :: k_na20__ne20__weak__wc12   = 158
  integer, parameter :: k_na20__p_ne19   = 159
  integer, parameter :: k_na20__he4_f16   = 160
  integer, parameter :: k_na20__he4_o16__weak__wc12   = 161
  integer, parameter :: k_p_na20__mg21   = 162
  integer, parameter :: k_he4_na20__al24   = 163
  integer, parameter :: k_p_na20__he4_ne17   = 164
  integer, parameter :: k_he4_na20__p_mg23   = 165
  integer, parameter :: k_mg23__na23__weak__wc12   = 166
  integer, parameter :: k_mg23__p_na22   = 167
  integer, parameter :: k_mg23__he4_ne19   = 168
  integer, parameter :: k_p_mg23__al24   = 169
  integer, parameter :: k_he4_mg23__si27   = 170
  integer, parameter :: k_p_mg23__he4_na20   = 171
  integer, parameter :: k_he4_mg23__p_al26   = 172
  integer, parameter :: k_mg24__p_na23   = 173
  integer, parameter :: k_mg24__he4_ne20   = 174
  integer, parameter :: k_p_mg24__al25   = 175
  integer, parameter :: k_he4_mg24__si28   = 176
  integer, parameter :: k_p_mg24__he4_na21   = 177
  integer, parameter :: k_he4_mg24__p_al27   = 178
  integer, parameter :: k_mg25__p_na24   = 179
  integer, parameter :: k_mg25__he4_ne21   = 180
  integer, parameter :: k_p_mg25__al26   = 181
  integer, parameter :: k_he4_mg25__si29   = 182
  integer, parameter :: k_p_mg25__he4_na22   = 183
  integer, parameter :: k_he4_mg25__p_al28   = 184
  integer, parameter :: k_mg26__he4_ne22   = 185
  integer, parameter :: k_p_mg26__al27   = 186
  integer, parameter :: k_he4_mg26__si30   = 187
  integer, parameter :: k_p_mg26__he4_na23   = 188
  integer, parameter :: k_al25__mg25__weak__wc12   = 189
  integer, parameter :: k_al25__p_mg24   = 190
  integer, parameter :: k_al25__he4_na21   = 191
  integer, parameter :: k_p_al25__si26   = 192
  integer, parameter :: k_he4_al25__p29   = 193
  integer, parameter :: k_p_al25__he4_mg22   = 194
  integer, parameter :: k_he4_al25__p_si28   = 195
  integer, parameter :: k_al26__mg26__weak__wc12   = 196
  integer, parameter :: k_al26__p_mg25   = 197
  integer, parameter :: k_al26__he4_na22   = 198
  integer, parameter :: k_p_al26__si27   = 199
  integer, parameter :: k_he4_al26__p30   = 200
  integer, parameter :: k_p_al26__he4_mg23   = 201
  integer, parameter :: k_he4_al26__p_si29   = 202
  integer, parameter :: k_al27__p_mg26   = 203
  integer, parameter :: k_al27__he4_na23   = 204
  integer, parameter :: k_p_al27__si28   = 205
  integer, parameter :: k_he4_al27__p31   = 206
  integer, parameter :: k_p_al27__he4_mg24   = 207
  integer, parameter :: k_he4_al27__p_si30   = 208
  integer, parameter :: k_al28__he4_na24   = 209
  integer, parameter :: k_p_al28__si29   = 210
  integer, parameter :: k_he4_al28__p32   = 211
  integer, parameter :: k_p_al28__he4_mg25   = 212
  integer, parameter :: k_al23__mg23__weak__wc12   = 213
  integer, parameter :: k_al23__p_mg22   = 214
  integer, parameter :: k_al23__p_na22__weak__wc12   = 215
  integer, parameter :: k_al23__he4_na19   = 216
  integer, parameter :: k_p_al23__si24   = 217
  integer, parameter :: k_he4_al23__p27   = 218
  integer, parameter :: k_he4_al23__p_si26   = 219
  integer, parameter :: k_si26__al26__weak__wc12   = 220
  integer, parameter :: k_si26__p_al25   = 221
  integer, parameter :: k_si26__he4_mg22   = 222
  integer, parameter :: k_p_si26__p27   = 223
  integer, parameter :: k_he4_si26__s30   = 224
  integer, parameter :: k_p_si26__he4_al23   = 225
  integer, parameter :: k_he4_si26__p_p29   = 226
  integer, parameter :: k_ne17__f17__weak__wc17   = 227
  integer, parameter :: k_ne17__p_o16__weak__wc12   = 228
  integer, parameter :: k_he4_ne17__mg21   = 229
  integer, parameter :: k_he4_ne17__p_na20   = 230
  integer, parameter :: k_mg21__na21__weak__wc12   = 231
  integer, parameter :: k_mg21__p_ne20__weak__wc12   = 232
  integer, parameter :: k_mg21__p_na20   = 233
  integer, parameter :: k_mg21__he4_ne17   = 234
  integer, parameter :: k_p_mg21__al22   = 235
  integer, parameter :: k_he4_mg21__si25   = 236
  integer, parameter :: k_he4_mg21__p_al24   = 237
  integer, parameter :: k_al24__mg24__weak__wc12   = 238
  integer, parameter :: k_al24__p_mg23   = 239
  integer, parameter :: k_al24__p_na23__weak__wc12   = 240
  integer, parameter :: k_al24__he4_na20   = 241
  integer, parameter :: k_al24__he4_ne20__weak__wc12   = 242
  integer, parameter :: k_p_al24__si25   = 243
  integer, parameter :: k_he4_al24__p28   = 244
  integer, parameter :: k_p_al24__he4_mg21   = 245
  integer, parameter :: k_he4_al24__p_si27   = 246
  integer, parameter :: k_si27__al27__weak__wc12   = 247
  integer, parameter :: k_si27__p_al26   = 248
  integer, parameter :: k_si27__he4_mg23   = 249
  integer, parameter :: k_p_si27__p28   = 250
  integer, parameter :: k_he4_si27__s31   = 251
  integer, parameter :: k_p_si27__he4_al24   = 252
  integer, parameter :: k_he4_si27__p_p30   = 253
  integer, parameter :: k_si28__p_al27   = 254
  integer, parameter :: k_si28__he4_mg24   = 255
  integer, parameter :: k_p_si28__p29   = 256
  integer, parameter :: k_he4_si28__s32   = 257
  integer, parameter :: k_p_si28__he4_al25   = 258
  integer, parameter :: k_he4_si28__p_p31   = 259
  integer, parameter :: k_si29__p_al28   = 260
  integer, parameter :: k_si29__he4_mg25   = 261
  integer, parameter :: k_p_si29__p30   = 262
  integer, parameter :: k_he4_si29__s33   = 263
  integer, parameter :: k_p_si29__he4_al26   = 264
  integer, parameter :: k_he4_si29__p_p32   = 265
  integer, parameter :: k_si30__he4_mg26   = 266
  integer, parameter :: k_p_si30__p31   = 267
  integer, parameter :: k_he4_si30__s34   = 268
  integer, parameter :: k_p_si30__he4_al27   = 269
  integer, parameter :: k_he4_si30__p_p33   = 270
  integer, parameter :: k_p29__si29__weak__wc12   = 271
  integer, parameter :: k_p29__p_si28   = 272
  integer, parameter :: k_p29__he4_al25   = 273
  integer, parameter :: k_p_p29__s30   = 274
  integer, parameter :: k_he4_p29__cl33   = 275
  integer, parameter :: k_p_p29__he4_si26   = 276
  integer, parameter :: k_he4_p29__p_s32   = 277
  integer, parameter :: k_p30__si30__weak__wc12   = 278
  integer, parameter :: k_p30__p_si29   = 279
  integer, parameter :: k_p30__he4_al26   = 280
  integer, parameter :: k_p_p30__s31   = 281
  integer, parameter :: k_he4_p30__cl34   = 282
  integer, parameter :: k_p_p30__he4_si27   = 283
  integer, parameter :: k_he4_p30__p_s33   = 284
  integer, parameter :: k_p31__p_si30   = 285
  integer, parameter :: k_p31__he4_al27   = 286
  integer, parameter :: k_p_p31__s32   = 287
  integer, parameter :: k_he4_p31__cl35   = 288
  integer, parameter :: k_p_p31__he4_si28   = 289
  integer, parameter :: k_he4_p31__p_s34   = 290
  integer, parameter :: k_p32__he4_al28   = 291
  integer, parameter :: k_p_p32__s33   = 292
  integer, parameter :: k_he4_p32__cl36   = 293
  integer, parameter :: k_p_p32__he4_si29   = 294
  integer, parameter :: k_he4_p32__p_s35   = 295
  integer, parameter :: k_si24__al24__weak__wc12   = 296
  integer, parameter :: k_si24__p_al23   = 297
  integer, parameter :: k_si24__p_mg23__weak__wc12   = 298
  integer, parameter :: k_he4_si24__s28   = 299
  integer, parameter :: k_he4_si24__p_p27   = 300
  integer, parameter :: k_p27__si27__weak__wc12   = 301
  integer, parameter :: k_p27__p_si26   = 302
  integer, parameter :: k_p27__p_al26__weak__wc12   = 303
  integer, parameter :: k_p27__he4_al23   = 304
  integer, parameter :: k_p_p27__s28   = 305
  integer, parameter :: k_he4_p27__cl31   = 306
  integer, parameter :: k_p_p27__he4_si24   = 307
  integer, parameter :: k_he4_p27__p_s30   = 308
  integer, parameter :: k_s30__p30__weak__wc12   = 309
  integer, parameter :: k_s30__p_p29   = 310
  integer, parameter :: k_s30__he4_si26   = 311
  integer, parameter :: k_p_s30__cl31   = 312
  integer, parameter :: k_he4_s30__ar34   = 313
  integer, parameter :: k_p_s30__he4_p27   = 314
  integer, parameter :: k_he4_s30__p_cl33   = 315
  integer, parameter :: k_al22__mg22__weak__wc12   = 316
  integer, parameter :: k_al22__p_mg21   = 317
  integer, parameter :: k_al22__p_na21__weak__wc17   = 318
  integer, parameter :: k_al22__he4_ne18__weak__wc12   = 319
  integer, parameter :: k_al22__p_p_ne20__weak__wc12   = 320
  integer, parameter :: k_he4_al22__p26   = 321
  integer, parameter :: k_he4_al22__p_si25   = 322
  integer, parameter :: k_si25__al25__weak__wc12   = 323
  integer, parameter :: k_si25__p_al24   = 324
  integer, parameter :: k_si25__p_mg24__weak__wc12   = 325
  integer, parameter :: k_si25__he4_mg21   = 326
  integer, parameter :: k_p_si25__p26   = 327
  integer, parameter :: k_he4_si25__s29   = 328
  integer, parameter :: k_p_si25__he4_al22   = 329
  integer, parameter :: k_he4_si25__p_p28   = 330
  integer, parameter :: k_p28__si28__weak__wc12   = 331
  integer, parameter :: k_p28__p_si27   = 332
  integer, parameter :: k_p28__p_al27__weak__wc12   = 333
  integer, parameter :: k_p28__he4_al24   = 334
  integer, parameter :: k_p28__he4_mg24__weak__wc12   = 335
  integer, parameter :: k_p_p28__s29   = 336
  integer, parameter :: k_he4_p28__cl32   = 337
  integer, parameter :: k_p_p28__he4_si25   = 338
  integer, parameter :: k_he4_p28__p_s31   = 339
  integer, parameter :: k_s31__p31__weak__wc12   = 340
  integer, parameter :: k_s31__p_p30   = 341
  integer, parameter :: k_s31__he4_si27   = 342
  integer, parameter :: k_p_s31__cl32   = 343
  integer, parameter :: k_he4_s31__ar35   = 344
  integer, parameter :: k_p_s31__he4_p28   = 345
  integer, parameter :: k_he4_s31__p_cl34   = 346
  integer, parameter :: k_s32__p_p31   = 347
  integer, parameter :: k_s32__he4_si28   = 348
  integer, parameter :: k_p_s32__cl33   = 349
  integer, parameter :: k_he4_s32__ar36   = 350
  integer, parameter :: k_p_s32__he4_p29   = 351
  integer, parameter :: k_he4_s32__p_cl35   = 352
  integer, parameter :: k_s33__p_p32   = 353
  integer, parameter :: k_s33__he4_si29   = 354
  integer, parameter :: k_p_s33__cl34   = 355
  integer, parameter :: k_he4_s33__ar37   = 356
  integer, parameter :: k_p_s33__he4_p30   = 357
  integer, parameter :: k_he4_s33__p_cl36   = 358
  integer, parameter :: k_p_p33__s34   = 359
  integer, parameter :: k_he4_p33__cl37   = 360
  integer, parameter :: k_p_p33__he4_si30   = 361
  integer, parameter :: k_s34__p_p33   = 362
  integer, parameter :: k_s34__he4_si30   = 363
  integer, parameter :: k_p_s34__cl35   = 364
  integer, parameter :: k_he4_s34__ar38   = 365
  integer, parameter :: k_p_s34__he4_p31   = 366
  integer, parameter :: k_he4_s34__p_cl37   = 367
  integer, parameter :: k_cl33__s33__weak__wc12   = 368
  integer, parameter :: k_cl33__p_s32   = 369
  integer, parameter :: k_cl33__he4_p29   = 370
  integer, parameter :: k_p_cl33__ar34   = 371
  integer, parameter :: k_he4_cl33__k37   = 372
  integer, parameter :: k_p_cl33__he4_s30   = 373
  integer, parameter :: k_he4_cl33__p_ar36   = 374
  integer, parameter :: k_cl34__s34__weak__wc12   = 375
  integer, parameter :: k_cl34__p_s33   = 376
  integer, parameter :: k_cl34__he4_p30   = 377
  integer, parameter :: k_p_cl34__ar35   = 378
  integer, parameter :: k_he4_cl34__k38   = 379
  integer, parameter :: k_p_cl34__he4_s31   = 380
  integer, parameter :: k_he4_cl34__p_ar37   = 381
  integer, parameter :: k_cl35__p_s34   = 382
  integer, parameter :: k_cl35__he4_p31   = 383
  integer, parameter :: k_p_cl35__ar36   = 384
  integer, parameter :: k_he4_cl35__k39   = 385
  integer, parameter :: k_p_cl35__he4_s32   = 386
  integer, parameter :: k_he4_cl35__p_ar38   = 387
  integer, parameter :: k_p_s35__cl36   = 388
  integer, parameter :: k_he4_s35__ar39   = 389
  integer, parameter :: k_p_s35__he4_p32   = 390
  integer, parameter :: k_cl36__p_s35   = 391
  integer, parameter :: k_cl36__he4_p32   = 392
  integer, parameter :: k_p_cl36__ar37   = 393
  integer, parameter :: k_he4_cl36__k40   = 394
  integer, parameter :: k_p_cl36__he4_s33   = 395
  integer, parameter :: k_he4_cl36__p_ar39   = 396
  integer, parameter :: k_s28__p28__weak__wc12   = 397
  integer, parameter :: k_s28__p_p27   = 398
  integer, parameter :: k_s28__p_si27__weak__wc12   = 399
  integer, parameter :: k_s28__he4_si24   = 400
  integer, parameter :: k_p_s28__cl29   = 401
  integer, parameter :: k_he4_s28__ar32   = 402
  integer, parameter :: k_he4_s28__p_cl31   = 403
  integer, parameter :: k_cl31__s31__weak__wc12   = 404
  integer, parameter :: k_cl31__p_s30   = 405
  integer, parameter :: k_cl31__p_p30__weak__wc12   = 406
  integer, parameter :: k_cl31__he4_p27   = 407
  integer, parameter :: k_p_cl31__ar32   = 408
  integer, parameter :: k_he4_cl31__k35   = 409
  integer, parameter :: k_p_cl31__he4_s28   = 410
  integer, parameter :: k_he4_cl31__p_ar34   = 411
  integer, parameter :: k_ar34__cl34__weak__wc12   = 412
  integer, parameter :: k_ar34__p_cl33   = 413
  integer, parameter :: k_ar34__he4_s30   = 414
  integer, parameter :: k_p_ar34__k35   = 415
  integer, parameter :: k_he4_ar34__ca38   = 416
  integer, parameter :: k_p_ar34__he4_cl31   = 417
  integer, parameter :: k_he4_ar34__p_k37   = 418
  integer, parameter :: k_p26__si26__weak__wc12   = 419
  integer, parameter :: k_p26__p_si25   = 420
  integer, parameter :: k_p26__he4_al22   = 421
  integer, parameter :: k_he4_p26__cl30   = 422
  integer, parameter :: k_he4_p26__p_s29   = 423
  integer, parameter :: k_s29__p29__weak__wc12   = 424
  integer, parameter :: k_s29__p_p28   = 425
  integer, parameter :: k_s29__p_si28__weak__wc12   = 426
  integer, parameter :: k_s29__he4_si25   = 427
  integer, parameter :: k_p_s29__cl30   = 428
  integer, parameter :: k_he4_s29__ar33   = 429
  integer, parameter :: k_p_s29__he4_p26   = 430
  integer, parameter :: k_he4_s29__p_cl32   = 431
  integer, parameter :: k_cl32__s32__weak__wc12   = 432
  integer, parameter :: k_cl32__p_s31   = 433
  integer, parameter :: k_cl32__p_p31__weak__wc12   = 434
  integer, parameter :: k_cl32__he4_p28   = 435
  integer, parameter :: k_cl32__he4_si28__weak__wc12   = 436
  integer, parameter :: k_p_cl32__ar33   = 437
  integer, parameter :: k_he4_cl32__k36   = 438
  integer, parameter :: k_p_cl32__he4_s29   = 439
  integer, parameter :: k_he4_cl32__p_ar35   = 440
  integer, parameter :: k_ar35__cl35__weak__wc12   = 441
  integer, parameter :: k_ar35__p_cl34   = 442
  integer, parameter :: k_ar35__he4_s31   = 443
  integer, parameter :: k_p_ar35__k36   = 444
  integer, parameter :: k_he4_ar35__ca39   = 445
  integer, parameter :: k_p_ar35__he4_cl32   = 446
  integer, parameter :: k_he4_ar35__p_k38   = 447
  integer, parameter :: k_ar36__p_cl35   = 448
  integer, parameter :: k_ar36__he4_s32   = 449
  integer, parameter :: k_p_ar36__k37   = 450
  integer, parameter :: k_he4_ar36__ca40   = 451
  integer, parameter :: k_p_ar36__he4_cl33   = 452
  integer, parameter :: k_he4_ar36__p_k39   = 453
  integer, parameter :: k_ar37__cl37__weak__wc12   = 454
  integer, parameter :: k_ar37__p_cl36   = 455
  integer, parameter :: k_ar37__he4_s33   = 456
  integer, parameter :: k_p_ar37__k38   = 457
  integer, parameter :: k_he4_ar37__ca41   = 458
  integer, parameter :: k_p_ar37__he4_cl34   = 459
  integer, parameter :: k_he4_ar37__p_k40   = 460
  integer, parameter :: k_cl37__he4_p33   = 461
  integer, parameter :: k_p_cl37__ar38   = 462
  integer, parameter :: k_he4_cl37__k41   = 463
  integer, parameter :: k_p_cl37__he4_s34   = 464
  integer, parameter :: k_ar38__p_cl37   = 465
  integer, parameter :: k_ar38__he4_s34   = 466
  integer, parameter :: k_p_ar38__k39   = 467
  integer, parameter :: k_he4_ar38__ca42   = 468
  integer, parameter :: k_p_ar38__he4_cl35   = 469
  integer, parameter :: k_he4_ar38__p_k41   = 470
  integer, parameter :: k_k37__ar37__weak__wc12   = 471
  integer, parameter :: k_k37__p_ar36   = 472
  integer, parameter :: k_k37__he4_cl33   = 473
  integer, parameter :: k_p_k37__ca38   = 474
  integer, parameter :: k_he4_k37__sc41   = 475
  integer, parameter :: k_p_k37__he4_ar34   = 476
  integer, parameter :: k_he4_k37__p_ca40   = 477
  integer, parameter :: k_k38__ar38__weak__wc12   = 478
  integer, parameter :: k_k38__p_ar37   = 479
  integer, parameter :: k_k38__he4_cl34   = 480
  integer, parameter :: k_p_k38__ca39   = 481
  integer, parameter :: k_he4_k38__sc42   = 482
  integer, parameter :: k_p_k38__he4_ar35   = 483
  integer, parameter :: k_he4_k38__p_ca41   = 484
  integer, parameter :: k_k39__p_ar38   = 485
  integer, parameter :: k_k39__he4_cl35   = 486
  integer, parameter :: k_p_k39__ca40   = 487
  integer, parameter :: k_he4_k39__sc43   = 488
  integer, parameter :: k_p_k39__he4_ar36   = 489
  integer, parameter :: k_he4_k39__p_ca42   = 490
  integer, parameter :: k_ar39__he4_s35   = 491
  integer, parameter :: k_p_ar39__k40   = 492
  integer, parameter :: k_he4_ar39__ca43   = 493
  integer, parameter :: k_p_ar39__he4_cl36   = 494
  integer, parameter :: k_k40__p_ar39   = 495
  integer, parameter :: k_k40__he4_cl36   = 496
  integer, parameter :: k_p_k40__ca41   = 497
  integer, parameter :: k_he4_k40__sc44   = 498
  integer, parameter :: k_p_k40__he4_ar37   = 499
  integer, parameter :: k_he4_k40__p_ca43   = 500
  integer, parameter :: k_cl29__s29__weak__bqa_pos_   = 501
  integer, parameter :: k_cl29__p_s28   = 502
  integer, parameter :: k_he4_cl29__k33   = 503
  integer, parameter :: k_he4_cl29__p_ar32   = 504
  integer, parameter :: k_ar32__cl32__weak__wc12   = 505
  integer, parameter :: k_ar32__p_cl31   = 506
  integer, parameter :: k_ar32__p_s31__weak__wc12   = 507
  integer, parameter :: k_ar32__he4_s28   = 508
  integer, parameter :: k_p_ar32__k33   = 509
  integer, parameter :: k_he4_ar32__ca36   = 510
  integer, parameter :: k_p_ar32__he4_cl29   = 511
  integer, parameter :: k_he4_ar32__p_k35   = 512
  integer, parameter :: k_k35__ar35__weak__wc12   = 513
  integer, parameter :: k_k35__p_ar34   = 514
  integer, parameter :: k_k35__p_cl34__weak__wc12   = 515
  integer, parameter :: k_k35__he4_cl31   = 516
  integer, parameter :: k_p_k35__ca36   = 517
  integer, parameter :: k_he4_k35__sc39   = 518
  integer, parameter :: k_p_k35__he4_ar32   = 519
  integer, parameter :: k_he4_k35__p_ca38   = 520
  integer, parameter :: k_ca38__k38__weak__wc12   = 521
  integer, parameter :: k_ca38__p_k37   = 522
  integer, parameter :: k_ca38__he4_ar34   = 523
  integer, parameter :: k_p_ca38__sc39   = 524
  integer, parameter :: k_he4_ca38__ti42   = 525
  integer, parameter :: k_p_ca38__he4_k35   = 526
  integer, parameter :: k_he4_ca38__p_sc41   = 527
  integer, parameter :: k_p_p_ca38__he4_ca36   = 528
  integer, parameter :: k_cl30__s30__weak__bqa_pos_   = 529
  integer, parameter :: k_cl30__p_s29   = 530
  integer, parameter :: k_cl30__he4_p26   = 531
  integer, parameter :: k_p_cl30__ar31   = 532
  integer, parameter :: k_he4_cl30__k34   = 533
  integer, parameter :: k_he4_cl30__p_ar33   = 534
  integer, parameter :: k_ar33__cl33__weak__wc12   = 535
  integer, parameter :: k_ar33__p_cl32   = 536
  integer, parameter :: k_ar33__p_s32__weak__wc12   = 537
  integer, parameter :: k_ar33__he4_s29   = 538
  integer, parameter :: k_p_ar33__k34   = 539
  integer, parameter :: k_he4_ar33__ca37   = 540
  integer, parameter :: k_p_ar33__he4_cl30   = 541
  integer, parameter :: k_he4_ar33__p_k36   = 542
  integer, parameter :: k_k36__ar36__weak__wc12   = 543
  integer, parameter :: k_k36__p_ar35   = 544
  integer, parameter :: k_k36__p_cl35__weak__wc12   = 545
  integer, parameter :: k_k36__he4_cl32   = 546
  integer, parameter :: k_k36__he4_s32__weak__wc12   = 547
  integer, parameter :: k_p_k36__ca37   = 548
  integer, parameter :: k_he4_k36__sc40   = 549
  integer, parameter :: k_p_k36__he4_ar33   = 550
  integer, parameter :: k_he4_k36__p_ca39   = 551
  integer, parameter :: k_ca39__k39__weak__wc12   = 552
  integer, parameter :: k_ca39__p_k38   = 553
  integer, parameter :: k_ca39__he4_ar35   = 554
  integer, parameter :: k_p_ca39__sc40   = 555
  integer, parameter :: k_he4_ca39__ti43   = 556
  integer, parameter :: k_p_ca39__he4_k36   = 557
  integer, parameter :: k_he4_ca39__p_sc42   = 558
  integer, parameter :: k_ca40__p_k39   = 559
  integer, parameter :: k_ca40__he4_ar36   = 560
  integer, parameter :: k_p_ca40__sc41   = 561
  integer, parameter :: k_he4_ca40__ti44   = 562
  integer, parameter :: k_p_ca40__he4_k37   = 563
  integer, parameter :: k_he4_ca40__p_sc43   = 564
  integer, parameter :: k_ca41__k41__weak__wc12   = 565
  integer, parameter :: k_ca41__p_k40   = 566
  integer, parameter :: k_ca41__he4_ar37   = 567
  integer, parameter :: k_p_ca41__sc42   = 568
  integer, parameter :: k_he4_ca41__ti45   = 569
  integer, parameter :: k_p_ca41__he4_k38   = 570
  integer, parameter :: k_he4_ca41__p_sc44   = 571
  integer, parameter :: k_k41__he4_cl37   = 572
  integer, parameter :: k_p_k41__ca42   = 573
  integer, parameter :: k_he4_k41__sc45   = 574
  integer, parameter :: k_p_k41__he4_ar38   = 575
  integer, parameter :: k_he4_k41__p_ca44   = 576
  integer, parameter :: k_ca42__p_k41   = 577
  integer, parameter :: k_ca42__he4_ar38   = 578
  integer, parameter :: k_p_ca42__sc43   = 579
  integer, parameter :: k_he4_ca42__ti46   = 580
  integer, parameter :: k_p_ca42__he4_k39   = 581
  integer, parameter :: k_he4_ca42__p_sc45   = 582
  integer, parameter :: k_sc41__ca41__weak__wc12   = 583
  integer, parameter :: k_sc41__p_ca40   = 584
  integer, parameter :: k_sc41__he4_k37   = 585
  integer, parameter :: k_p_sc41__ti42   = 586
  integer, parameter :: k_he4_sc41__v45   = 587
  integer, parameter :: k_p_sc41__he4_ca38   = 588
  integer, parameter :: k_he4_sc41__p_ti44   = 589
  integer, parameter :: k_sc42__ca42__weak__wc12   = 590
  integer, parameter :: k_sc42__p_ca41   = 591
  integer, parameter :: k_sc42__he4_k38   = 592
  integer, parameter :: k_p_sc42__ti43   = 593
  integer, parameter :: k_he4_sc42__v46   = 594
  integer, parameter :: k_p_sc42__he4_ca39   = 595
  integer, parameter :: k_he4_sc42__p_ti45   = 596
  integer, parameter :: k_sc43__ca43__weak__wc12   = 597
  integer, parameter :: k_sc43__p_ca42   = 598
  integer, parameter :: k_sc43__he4_k39   = 599
  integer, parameter :: k_p_sc43__ti44   = 600
  integer, parameter :: k_he4_sc43__v47   = 601
  integer, parameter :: k_p_sc43__he4_ca40   = 602
  integer, parameter :: k_he4_sc43__p_ti46   = 603
  integer, parameter :: k_ca43__he4_ar39   = 604
  integer, parameter :: k_p_ca43__sc44   = 605
  integer, parameter :: k_he4_ca43__ti47   = 606
  integer, parameter :: k_p_ca43__he4_k40   = 607
  integer, parameter :: k_he4_ca43__p_sc46   = 608
  integer, parameter :: k_sc44__ca44__weak__wc12   = 609
  integer, parameter :: k_sc44__p_ca43   = 610
  integer, parameter :: k_sc44__he4_k40   = 611
  integer, parameter :: k_p_sc44__ti45   = 612
  integer, parameter :: k_he4_sc44__v48   = 613
  integer, parameter :: k_p_sc44__he4_ca41   = 614
  integer, parameter :: k_he4_sc44__p_ti47   = 615
  integer, parameter :: k_k33__ar33__weak__bqa_pos_   = 616
  integer, parameter :: k_k33__p_ar32   = 617
  integer, parameter :: k_k33__he4_cl29   = 618
  integer, parameter :: k_p_k33__ca34   = 619
  integer, parameter :: k_he4_k33__sc37   = 620
  integer, parameter :: k_he4_k33__p_ca36   = 621
  integer, parameter :: k_ca36__k36__weak__wc12   = 622
  integer, parameter :: k_ca36__p_k35   = 623
  integer, parameter :: k_ca36__p_ar35__weak__wc12   = 624
  integer, parameter :: k_ca36__he4_ar32   = 625
  integer, parameter :: k_p_ca36__sc37   = 626
  integer, parameter :: k_he4_ca36__ti40   = 627
  integer, parameter :: k_p_ca36__he4_k33   = 628
  integer, parameter :: k_he4_ca36__p_sc39   = 629
  integer, parameter :: k_he4_ca36__p_p_ca38   = 630
  integer, parameter :: k_sc39__ca39__weak__mo97   = 631
  integer, parameter :: k_sc39__p_ca38   = 632
  integer, parameter :: k_sc39__he4_k35   = 633
  integer, parameter :: k_p_sc39__ti40   = 634
  integer, parameter :: k_he4_sc39__v43   = 635
  integer, parameter :: k_p_sc39__he4_ca36   = 636
  integer, parameter :: k_he4_sc39__p_ti42   = 637
  integer, parameter :: k_ti42__sc42__weak__wc12   = 638
  integer, parameter :: k_ti42__p_sc41   = 639
  integer, parameter :: k_ti42__he4_ca38   = 640
  integer, parameter :: k_p_ti42__v43   = 641
  integer, parameter :: k_he4_ti42__cr46   = 642
  integer, parameter :: k_p_ti42__he4_sc39   = 643
  integer, parameter :: k_he4_ti42__p_v45   = 644
  integer, parameter :: k_ar31__cl31__weak__wc12   = 645
  integer, parameter :: k_ar31__p_cl30   = 646
  integer, parameter :: k_ar31__p_s30__weak__wc17   = 647
  integer, parameter :: k_ar31__p_p_p29__weak__wc12   = 648
  integer, parameter :: k_he4_ar31__ca35   = 649
  integer, parameter :: k_he4_ar31__p_k34   = 650
  integer, parameter :: k_k34__ar34__weak__bqa_pos_   = 651
  integer, parameter :: k_k34__p_ar33   = 652
  integer, parameter :: k_k34__he4_cl30   = 653
  integer, parameter :: k_p_k34__ca35   = 654
  integer, parameter :: k_he4_k34__sc38   = 655
  integer, parameter :: k_p_k34__he4_ar31   = 656
  integer, parameter :: k_he4_k34__p_ca37   = 657
  integer, parameter :: k_ca37__k37__weak__wc12   = 658
  integer, parameter :: k_ca37__p_ar36__weak__wc12   = 659
  integer, parameter :: k_ca37__p_k36   = 660
  integer, parameter :: k_ca37__he4_ar33   = 661
  integer, parameter :: k_p_ca37__sc38   = 662
  integer, parameter :: k_he4_ca37__ti41   = 663
  integer, parameter :: k_p_ca37__he4_k34   = 664
  integer, parameter :: k_he4_ca37__p_sc40   = 665
  integer, parameter :: k_sc40__ca40__weak__wc12   = 666
  integer, parameter :: k_sc40__p_ca39   = 667
  integer, parameter :: k_sc40__p_k39__weak__wc12   = 668
  integer, parameter :: k_sc40__he4_k36   = 669
  integer, parameter :: k_sc40__he4_ar36__weak__wc12   = 670
  integer, parameter :: k_p_sc40__ti41   = 671
  integer, parameter :: k_he4_sc40__v44   = 672
  integer, parameter :: k_p_sc40__he4_ca37   = 673
  integer, parameter :: k_he4_sc40__p_ti43   = 674
  integer, parameter :: k_ti43__sc43__weak__wc12   = 675
  integer, parameter :: k_ti43__p_sc42   = 676
  integer, parameter :: k_ti43__he4_ca39   = 677
  integer, parameter :: k_p_ti43__v44   = 678
  integer, parameter :: k_he4_ti43__cr47   = 679
  integer, parameter :: k_p_ti43__he4_sc40   = 680
  integer, parameter :: k_he4_ti43__p_v46   = 681
  integer, parameter :: k_ti44__sc44__weak__wc12   = 682
  integer, parameter :: k_ti44__p_sc43   = 683
  integer, parameter :: k_ti44__he4_ca40   = 684
  integer, parameter :: k_p_ti44__v45   = 685
  integer, parameter :: k_he4_ti44__cr48   = 686
  integer, parameter :: k_p_ti44__he4_sc41   = 687
  integer, parameter :: k_he4_ti44__p_v47   = 688
  integer, parameter :: k_ti45__sc45__weak__wc12   = 689
  integer, parameter :: k_ti45__p_sc44   = 690
  integer, parameter :: k_ti45__he4_ca41   = 691
  integer, parameter :: k_p_ti45__v46   = 692
  integer, parameter :: k_he4_ti45__cr49   = 693
  integer, parameter :: k_p_ti45__he4_sc42   = 694
  integer, parameter :: k_he4_ti45__p_v48   = 695
  integer, parameter :: k_p_ca44__sc45   = 696
  integer, parameter :: k_he4_ca44__ti48   = 697
  integer, parameter :: k_p_ca44__he4_k41   = 698
  integer, parameter :: k_sc45__p_ca44   = 699
  integer, parameter :: k_sc45__he4_k41   = 700
  integer, parameter :: k_p_sc45__ti46   = 701
  integer, parameter :: k_he4_sc45__v49   = 702
  integer, parameter :: k_p_sc45__he4_ca42   = 703
  integer, parameter :: k_he4_sc45__p_ti48   = 704
  integer, parameter :: k_ti46__p_sc45   = 705
  integer, parameter :: k_ti46__he4_ca42   = 706
  integer, parameter :: k_p_ti46__v47   = 707
  integer, parameter :: k_he4_ti46__cr50   = 708
  integer, parameter :: k_p_ti46__he4_sc43   = 709
  integer, parameter :: k_he4_ti46__p_v49   = 710
  integer, parameter :: k_v45__ti45__weak__wc12   = 711
  integer, parameter :: k_v45__p_ti44   = 712
  integer, parameter :: k_v45__he4_sc41   = 713
  integer, parameter :: k_p_v45__cr46   = 714
  integer, parameter :: k_he4_v45__mn49   = 715
  integer, parameter :: k_p_v45__he4_ti42   = 716
  integer, parameter :: k_he4_v45__p_cr48   = 717
  integer, parameter :: k_v46__ti46__weak__wc12   = 718
  integer, parameter :: k_v46__p_ti45   = 719
  integer, parameter :: k_v46__he4_sc42   = 720
  integer, parameter :: k_p_v46__cr47   = 721
  integer, parameter :: k_he4_v46__mn50   = 722
  integer, parameter :: k_p_v46__he4_ti43   = 723
  integer, parameter :: k_he4_v46__p_cr49   = 724
  integer, parameter :: k_v47__ti47__weak__wc12   = 725
  integer, parameter :: k_v47__p_ti46   = 726
  integer, parameter :: k_v47__he4_sc43   = 727
  integer, parameter :: k_p_v47__cr48   = 728
  integer, parameter :: k_he4_v47__mn51   = 729
  integer, parameter :: k_p_v47__he4_ti44   = 730
  integer, parameter :: k_he4_v47__p_cr50   = 731
  integer, parameter :: k_p_sc46__ti47   = 732
  integer, parameter :: k_he4_sc46__v50   = 733
  integer, parameter :: k_p_sc46__he4_ca43   = 734
  integer, parameter :: k_ti47__p_sc46   = 735
  integer, parameter :: k_ti47__he4_ca43   = 736
  integer, parameter :: k_p_ti47__v48   = 737
  integer, parameter :: k_he4_ti47__cr51   = 738
  integer, parameter :: k_p_ti47__he4_sc44   = 739
  integer, parameter :: k_he4_ti47__p_v50   = 740
  integer, parameter :: k_v48__ti48__weak__wc12   = 741
  integer, parameter :: k_v48__p_ti47   = 742
  integer, parameter :: k_v48__he4_sc44   = 743
  integer, parameter :: k_p_v48__cr49   = 744
  integer, parameter :: k_he4_v48__mn52   = 745
  integer, parameter :: k_p_v48__he4_ti45   = 746
  integer, parameter :: k_he4_v48__p_cr51   = 747
  integer, parameter :: k_ca34__k34__weak__bqa_pos_   = 748
  integer, parameter :: k_ca34__p_k33   = 749
  integer, parameter :: k_he4_ca34__ti38   = 750
  integer, parameter :: k_he4_ca34__p_sc37   = 751
  integer, parameter :: k_sc37__ca37__weak__bqa_pos_   = 752
  integer, parameter :: k_sc37__p_ca36   = 753
  integer, parameter :: k_sc37__he4_k33   = 754
  integer, parameter :: k_p_sc37__ti38   = 755
  integer, parameter :: k_he4_sc37__v41   = 756
  integer, parameter :: k_p_sc37__he4_ca34   = 757
  integer, parameter :: k_he4_sc37__p_ti40   = 758
  integer, parameter :: k_ti40__sc40__weak__wc17   = 759
  integer, parameter :: k_ti40__p_sc39   = 760
  integer, parameter :: k_ti40__p_ca39__weak__wc12   = 761
  integer, parameter :: k_ti40__he4_ca36   = 762
  integer, parameter :: k_p_ti40__v41   = 763
  integer, parameter :: k_he4_ti40__cr44   = 764
  integer, parameter :: k_p_ti40__he4_sc37   = 765
  integer, parameter :: k_he4_ti40__p_v43   = 766
  integer, parameter :: k_v43__ti43__weak__wc12   = 767
  integer, parameter :: k_v43__p_ti42   = 768
  integer, parameter :: k_v43__he4_sc39   = 769
  integer, parameter :: k_p_v43__cr44   = 770
  integer, parameter :: k_he4_v43__mn47   = 771
  integer, parameter :: k_p_v43__he4_ti40   = 772
  integer, parameter :: k_he4_v43__p_cr46   = 773
  integer, parameter :: k_cr46__v46__weak__wc12   = 774
  integer, parameter :: k_cr46__p_v45   = 775
  integer, parameter :: k_cr46__he4_ti42   = 776
  integer, parameter :: k_p_cr46__mn47   = 777
  integer, parameter :: k_he4_cr46__fe50   = 778
  integer, parameter :: k_p_cr46__he4_v43   = 779
  integer, parameter :: k_he4_cr46__p_mn49   = 780
  integer, parameter :: k_ca35__k35__weak__wc17   = 781
  integer, parameter :: k_ca35__p_k34   = 782
  integer, parameter :: k_ca35__p_ar34__weak__wc17   = 783
  integer, parameter :: k_ca35__he4_ar31   = 784
  integer, parameter :: k_ca35__p_p_cl33__weak__wc12   = 785
  integer, parameter :: k_p_ca35__sc36   = 786
  integer, parameter :: k_he4_ca35__ti39   = 787
  integer, parameter :: k_he4_ca35__p_sc38   = 788
  integer, parameter :: k_sc38__ca38__weak__mo97   = 789
  integer, parameter :: k_sc38__p_ca37   = 790
  integer, parameter :: k_sc38__he4_k34   = 791
  integer, parameter :: k_p_sc38__ti39   = 792
  integer, parameter :: k_he4_sc38__v42   = 793
  integer, parameter :: k_p_sc38__he4_ca35   = 794
  integer, parameter :: k_he4_sc38__p_ti41   = 795
  integer, parameter :: k_ti41__sc41__weak__wc17   = 796
  integer, parameter :: k_ti41__p_sc40   = 797
  integer, parameter :: k_ti41__p_ca40__weak__wc12   = 798
  integer, parameter :: k_ti41__he4_ca37   = 799
  integer, parameter :: k_p_ti41__v42   = 800
  integer, parameter :: k_he4_ti41__cr45   = 801
  integer, parameter :: k_p_ti41__he4_sc38   = 802
  integer, parameter :: k_he4_ti41__p_v44   = 803
  integer, parameter :: k_v44__ti44__weak__wc12   = 804
  integer, parameter :: k_v44__p_ti43   = 805
  integer, parameter :: k_v44__he4_sc40   = 806
  integer, parameter :: k_p_v44__cr45   = 807
  integer, parameter :: k_he4_v44__mn48   = 808
  integer, parameter :: k_p_v44__he4_ti41   = 809
  integer, parameter :: k_he4_v44__p_cr47   = 810
  integer, parameter :: k_cr47__v47__weak__wc12   = 811
  integer, parameter :: k_cr47__p_v46   = 812
  integer, parameter :: k_cr47__he4_ti43   = 813
  integer, parameter :: k_p_cr47__mn48   = 814
  integer, parameter :: k_he4_cr47__fe51   = 815
  integer, parameter :: k_p_cr47__he4_v44   = 816
  integer, parameter :: k_he4_cr47__p_mn50   = 817
  integer, parameter :: k_cr48__v48__weak__wc12   = 818
  integer, parameter :: k_cr48__p_v47   = 819
  integer, parameter :: k_cr48__he4_ti44   = 820
  integer, parameter :: k_p_cr48__mn49   = 821
  integer, parameter :: k_he4_cr48__fe52   = 822
  integer, parameter :: k_p_cr48__he4_v45   = 823
  integer, parameter :: k_he4_cr48__p_mn51   = 824
  integer, parameter :: k_cr49__v49__weak__wc12   = 825
  integer, parameter :: k_cr49__p_v48   = 826
  integer, parameter :: k_cr49__he4_ti45   = 827
  integer, parameter :: k_p_cr49__mn50   = 828
  integer, parameter :: k_he4_cr49__fe53   = 829
  integer, parameter :: k_p_cr49__he4_v46   = 830
  integer, parameter :: k_he4_cr49__p_mn52   = 831
  integer, parameter :: k_ti48__he4_ca44   = 832
  integer, parameter :: k_p_ti48__v49   = 833
  integer, parameter :: k_he4_ti48__cr52   = 834
  integer, parameter :: k_p_ti48__he4_sc45   = 835
  integer, parameter :: k_v49__p_ti48   = 836
  integer, parameter :: k_v49__he4_sc45   = 837
  integer, parameter :: k_p_v49__cr50   = 838
  integer, parameter :: k_he4_v49__mn53   = 839
  integer, parameter :: k_p_v49__he4_ti46   = 840
  integer, parameter :: k_he4_v49__p_cr52   = 841
  integer, parameter :: k_cr50__p_v49   = 842
  integer, parameter :: k_cr50__he4_ti46   = 843
  integer, parameter :: k_p_cr50__mn51   = 844
  integer, parameter :: k_he4_cr50__fe54   = 845
  integer, parameter :: k_p_cr50__he4_v47   = 846
  integer, parameter :: k_he4_cr50__p_mn53   = 847
  integer, parameter :: k_mn49__cr49__weak__wc12   = 848
  integer, parameter :: k_mn49__p_cr48   = 849
  integer, parameter :: k_mn49__he4_v45   = 850
  integer, parameter :: k_p_mn49__fe50   = 851
  integer, parameter :: k_he4_mn49__co53   = 852
  integer, parameter :: k_p_mn49__he4_cr46   = 853
  integer, parameter :: k_he4_mn49__p_fe52   = 854
  integer, parameter :: k_mn50__cr50__weak__wc12   = 855
  integer, parameter :: k_mn50__p_cr49   = 856
  integer, parameter :: k_mn50__he4_v46   = 857
  integer, parameter :: k_p_mn50__fe51   = 858
  integer, parameter :: k_he4_mn50__co54   = 859
  integer, parameter :: k_p_mn50__he4_cr47   = 860
  integer, parameter :: k_he4_mn50__p_fe53   = 861
  integer, parameter :: k_mn51__cr51__weak__wc12   = 862
  integer, parameter :: k_mn51__p_cr50   = 863
  integer, parameter :: k_mn51__he4_v47   = 864
  integer, parameter :: k_p_mn51__fe52   = 865
  integer, parameter :: k_he4_mn51__co55   = 866
  integer, parameter :: k_p_mn51__he4_cr48   = 867
  integer, parameter :: k_he4_mn51__p_fe54   = 868
  integer, parameter :: k_v50__he4_sc46   = 869
  integer, parameter :: k_p_v50__cr51   = 870
  integer, parameter :: k_he4_v50__mn54   = 871
  integer, parameter :: k_p_v50__he4_ti47   = 872
  integer, parameter :: k_cr51__p_v50   = 873
  integer, parameter :: k_cr51__he4_ti47   = 874
  integer, parameter :: k_p_cr51__mn52   = 875
  integer, parameter :: k_he4_cr51__fe55   = 876
  integer, parameter :: k_p_cr51__he4_v48   = 877
  integer, parameter :: k_he4_cr51__p_mn54   = 878
  integer, parameter :: k_mn52__cr52__weak__wc12   = 879
  integer, parameter :: k_mn52__p_cr51   = 880
  integer, parameter :: k_mn52__he4_v48   = 881
  integer, parameter :: k_p_mn52__fe53   = 882
  integer, parameter :: k_he4_mn52__co56   = 883
  integer, parameter :: k_p_mn52__he4_cr49   = 884
  integer, parameter :: k_he4_mn52__p_fe55   = 885
  integer, parameter :: k_ti38__sc38__weak__mo97   = 886
  integer, parameter :: k_ti38__p_sc37   = 887
  integer, parameter :: k_ti38__he4_ca34   = 888
  integer, parameter :: k_he4_ti38__cr42   = 889
  integer, parameter :: k_he4_ti38__p_v41   = 890
  integer, parameter :: k_v41__ti41__weak__bqa_pos_   = 891
  integer, parameter :: k_v41__p_ti40   = 892
  integer, parameter :: k_v41__he4_sc37   = 893
  integer, parameter :: k_p_v41__cr42   = 894
  integer, parameter :: k_he4_v41__mn45   = 895
  integer, parameter :: k_p_v41__he4_ti38   = 896
  integer, parameter :: k_he4_v41__p_cr44   = 897
  integer, parameter :: k_cr44__v44__weak__wc12   = 898
  integer, parameter :: k_cr44__p_v43   = 899
  integer, parameter :: k_cr44__p_ti43__weak__wc12   = 900
  integer, parameter :: k_cr44__he4_ti40   = 901
  integer, parameter :: k_p_cr44__mn45   = 902
  integer, parameter :: k_he4_cr44__fe48   = 903
  integer, parameter :: k_p_cr44__he4_v41   = 904
  integer, parameter :: k_he4_cr44__p_mn47   = 905
  integer, parameter :: k_mn47__cr47__weak__wc17   = 906
  integer, parameter :: k_mn47__p_cr46   = 907
  integer, parameter :: k_mn47__he4_v43   = 908
  integer, parameter :: k_p_mn47__fe48   = 909
  integer, parameter :: k_he4_mn47__co51   = 910
  integer, parameter :: k_p_mn47__he4_cr44   = 911
  integer, parameter :: k_he4_mn47__p_fe50   = 912
  integer, parameter :: k_fe50__mn50__weak__wc12   = 913
  integer, parameter :: k_fe50__p_mn49   = 914
  integer, parameter :: k_fe50__he4_cr46   = 915
  integer, parameter :: k_p_fe50__co51   = 916
  integer, parameter :: k_he4_fe50__ni54   = 917
  integer, parameter :: k_p_fe50__he4_mn47   = 918
  integer, parameter :: k_he4_fe50__p_co53   = 919
  integer, parameter :: k_sc36__ca36__weak__bqa_pos_   = 920
  integer, parameter :: k_sc36__p_ca35   = 921
  integer, parameter :: k_he4_sc36__v40   = 922
  integer, parameter :: k_he4_sc36__p_ti39   = 923
  integer, parameter :: k_ti39__sc39__weak__wc17   = 924
  integer, parameter :: k_ti39__p_sc38   = 925
  integer, parameter :: k_ti39__p_ca38__weak__wc12   = 926
  integer, parameter :: k_ti39__he4_ca35   = 927
  integer, parameter :: k_p_ti39__v40   = 928
  integer, parameter :: k_he4_ti39__cr43   = 929
  integer, parameter :: k_p_ti39__he4_sc36   = 930
  integer, parameter :: k_he4_ti39__p_v42   = 931
  integer, parameter :: k_v42__ti42__weak__mo97   = 932
  integer, parameter :: k_v42__p_ti41   = 933
  integer, parameter :: k_v42__he4_sc38   = 934
  integer, parameter :: k_p_v42__cr43   = 935
  integer, parameter :: k_he4_v42__mn46   = 936
  integer, parameter :: k_p_v42__he4_ti39   = 937
  integer, parameter :: k_he4_v42__p_cr45   = 938
  integer, parameter :: k_cr45__v45__weak__wc12   = 939
  integer, parameter :: k_cr45__p_v44   = 940
  integer, parameter :: k_cr45__p_ti44__weak__wc12   = 941
  integer, parameter :: k_cr45__he4_ti41   = 942
  integer, parameter :: k_p_cr45__mn46   = 943
  integer, parameter :: k_he4_cr45__fe49   = 944
  integer, parameter :: k_p_cr45__he4_v42   = 945
  integer, parameter :: k_he4_cr45__p_mn48   = 946
  integer, parameter :: k_mn48__cr48__weak__wc12   = 947
  integer, parameter :: k_mn48__p_cr47   = 948
  integer, parameter :: k_mn48__p_v47__weak__wc12   = 949
  integer, parameter :: k_mn48__he4_v44   = 950
  integer, parameter :: k_p_mn48__fe49   = 951
  integer, parameter :: k_he4_mn48__co52   = 952
  integer, parameter :: k_p_mn48__he4_cr45   = 953
  integer, parameter :: k_he4_mn48__p_fe51   = 954
  integer, parameter :: k_fe51__mn51__weak__wc12   = 955
  integer, parameter :: k_fe51__p_mn50   = 956
  integer, parameter :: k_fe51__he4_cr47   = 957
  integer, parameter :: k_p_fe51__co52   = 958
  integer, parameter :: k_he4_fe51__ni55   = 959
  integer, parameter :: k_p_fe51__he4_mn48   = 960
  integer, parameter :: k_he4_fe51__p_co54   = 961
  integer, parameter :: k_fe52__mn52__weak__wc12   = 962
  integer, parameter :: k_fe52__p_mn51   = 963
  integer, parameter :: k_fe52__he4_cr48   = 964
  integer, parameter :: k_p_fe52__co53   = 965
  integer, parameter :: k_he4_fe52__ni56   = 966
  integer, parameter :: k_p_fe52__he4_mn49   = 967
  integer, parameter :: k_he4_fe52__p_co55   = 968
  integer, parameter :: k_fe53__mn53__weak__wc12   = 969
  integer, parameter :: k_fe53__p_mn52   = 970
  integer, parameter :: k_fe53__he4_cr49   = 971
  integer, parameter :: k_p_fe53__co54   = 972
  integer, parameter :: k_p_fe53__he4_mn50   = 973
  integer, parameter :: k_he4_fe53__p_co56   = 974
  integer, parameter :: k_cr52__he4_ti48   = 975
  integer, parameter :: k_p_cr52__mn53   = 976
  integer, parameter :: k_he4_cr52__fe56   = 977
  integer, parameter :: k_p_cr52__he4_v49   = 978
  integer, parameter :: k_he4_cr52__p_mn55   = 979
  integer, parameter :: k_mn53__p_cr52   = 980
  integer, parameter :: k_mn53__he4_v49   = 981
  integer, parameter :: k_p_mn53__fe54   = 982
  integer, parameter :: k_p_mn53__he4_cr50   = 983
  integer, parameter :: k_he4_mn53__p_fe56   = 984
  integer, parameter :: k_fe54__p_mn53   = 985
  integer, parameter :: k_fe54__he4_cr50   = 986
  integer, parameter :: k_p_fe54__co55   = 987
  integer, parameter :: k_p_fe54__he4_mn51   = 988
  integer, parameter :: k_co53__fe53__weak__wc12   = 989
  integer, parameter :: k_co53__p_fe52   = 990
  integer, parameter :: k_co53__he4_mn49   = 991
  integer, parameter :: k_p_co53__ni54   = 992
  integer, parameter :: k_p_co53__he4_fe50   = 993
  integer, parameter :: k_he4_co53__p_ni56   = 994
  integer, parameter :: k_co54__fe54__weak__wc12   = 995
  integer, parameter :: k_co54__p_fe53   = 996
  integer, parameter :: k_co54__he4_mn50   = 997
  integer, parameter :: k_p_co54__ni55   = 998
  integer, parameter :: k_p_co54__he4_fe51   = 999
  integer, parameter :: k_co55__fe55__weak__wc12   = 1000
  integer, parameter :: k_co55__p_fe54   = 1001
  integer, parameter :: k_co55__he4_mn51   = 1002
  integer, parameter :: k_p_co55__ni56   = 1003
  integer, parameter :: k_p_co55__he4_fe52   = 1004
  integer, parameter :: k_mn54__he4_v50   = 1005
  integer, parameter :: k_p_mn54__fe55   = 1006
  integer, parameter :: k_p_mn54__he4_cr51   = 1007
  integer, parameter :: k_fe55__mn55__weak__wc12   = 1008
  integer, parameter :: k_fe55__p_mn54   = 1009
  integer, parameter :: k_fe55__he4_cr51   = 1010
  integer, parameter :: k_p_fe55__co56   = 1011
  integer, parameter :: k_p_fe55__he4_mn52   = 1012
  integer, parameter :: k_co56__fe56__weak__wc12   = 1013
  integer, parameter :: k_co56__p_fe55   = 1014
  integer, parameter :: k_co56__he4_mn52   = 1015
  integer, parameter :: k_p_co56__he4_fe53   = 1016
  integer, parameter :: k_cr42__v42__weak__wc12   = 1017
  integer, parameter :: k_cr42__p_v41   = 1018
  integer, parameter :: k_cr42__p_ti41__weak__wc12   = 1019
  integer, parameter :: k_cr42__he4_ti38   = 1020
  integer, parameter :: k_he4_cr42__fe46   = 1021
  integer, parameter :: k_he4_cr42__p_mn45   = 1022
  integer, parameter :: k_mn45__cr45__weak__bqa_pos_   = 1023
  integer, parameter :: k_mn45__p_cr44   = 1024
  integer, parameter :: k_mn45__he4_v41   = 1025
  integer, parameter :: k_p_mn45__fe46   = 1026
  integer, parameter :: k_he4_mn45__co49   = 1027
  integer, parameter :: k_p_mn45__he4_cr42   = 1028
  integer, parameter :: k_he4_mn45__p_fe48   = 1029
  integer, parameter :: k_fe48__mn48__weak__wc12   = 1030
  integer, parameter :: k_fe48__p_mn47   = 1031
  integer, parameter :: k_fe48__p_cr47__weak__wc12   = 1032
  integer, parameter :: k_fe48__he4_cr44   = 1033
  integer, parameter :: k_p_fe48__co49   = 1034
  integer, parameter :: k_he4_fe48__ni52   = 1035
  integer, parameter :: k_p_fe48__he4_mn45   = 1036
  integer, parameter :: k_he4_fe48__p_co51   = 1037
  integer, parameter :: k_co51__fe51__weak__mo97   = 1038
  integer, parameter :: k_co51__p_fe50   = 1039
  integer, parameter :: k_co51__he4_mn47   = 1040
  integer, parameter :: k_p_co51__ni52   = 1041
  integer, parameter :: k_p_co51__he4_fe48   = 1042
  integer, parameter :: k_he4_co51__p_ni54   = 1043
  integer, parameter :: k_ni54__co54__weak__wc12   = 1044
  integer, parameter :: k_ni54__p_co53   = 1045
  integer, parameter :: k_ni54__he4_fe50   = 1046
  integer, parameter :: k_p_ni54__he4_co51   = 1047
  integer, parameter :: k_v40__ti40__weak__bqa_pos_   = 1048
  integer, parameter :: k_v40__p_ti39   = 1049
  integer, parameter :: k_v40__he4_sc36   = 1050
  integer, parameter :: k_he4_v40__mn44   = 1051
  integer, parameter :: k_he4_v40__p_cr43   = 1052
  integer, parameter :: k_cr43__v43__weak__wc12   = 1053
  integer, parameter :: k_cr43__p_v42   = 1054
  integer, parameter :: k_cr43__p_ti42__weak__wc17   = 1055
  integer, parameter :: k_cr43__he4_ti39   = 1056
  integer, parameter :: k_cr43__p_p_sc41__weak__wc12   = 1057
  integer, parameter :: k_p_cr43__mn44   = 1058
  integer, parameter :: k_he4_cr43__fe47   = 1059
  integer, parameter :: k_p_cr43__he4_v40   = 1060
  integer, parameter :: k_he4_cr43__p_mn46   = 1061
  integer, parameter :: k_cr43__p_p_p_ca40__weak__wc12   = 1062
  integer, parameter :: k_mn46__cr46__weak__wc12   = 1063
  integer, parameter :: k_mn46__p_cr45   = 1064
  integer, parameter :: k_mn46__p_v45__weak__wc12   = 1065
  integer, parameter :: k_mn46__he4_v42   = 1066
  integer, parameter :: k_p_mn46__fe47   = 1067
  integer, parameter :: k_he4_mn46__co50   = 1068
  integer, parameter :: k_p_mn46__he4_cr43   = 1069
  integer, parameter :: k_he4_mn46__p_fe49   = 1070
  integer, parameter :: k_fe49__mn49__weak__wc12   = 1071
  integer, parameter :: k_fe49__p_mn48   = 1072
  integer, parameter :: k_fe49__p_cr48__weak__wc12   = 1073
  integer, parameter :: k_fe49__he4_cr45   = 1074
  integer, parameter :: k_p_fe49__co50   = 1075
  integer, parameter :: k_he4_fe49__ni53   = 1076
  integer, parameter :: k_p_fe49__he4_mn46   = 1077
  integer, parameter :: k_he4_fe49__p_co52   = 1078
  integer, parameter :: k_co52__fe52__weak__wc12   = 1079
  integer, parameter :: k_co52__p_fe51   = 1080
  integer, parameter :: k_co52__he4_mn48   = 1081
  integer, parameter :: k_p_co52__ni53   = 1082
  integer, parameter :: k_p_co52__he4_fe49   = 1083
  integer, parameter :: k_he4_co52__p_ni55   = 1084
  integer, parameter :: k_ni55__co55__weak__wc12   = 1085
  integer, parameter :: k_ni55__p_co54   = 1086
  integer, parameter :: k_ni55__he4_fe51   = 1087
  integer, parameter :: k_p_ni55__he4_co52   = 1088
  integer, parameter :: k_ni56__co56__weak__wc12   = 1089
  integer, parameter :: k_ni56__p_co55   = 1090
  integer, parameter :: k_ni56__he4_fe52   = 1091
  integer, parameter :: k_p_ni56__he4_co53   = 1092
  integer, parameter :: k_p_mn55__fe56   = 1093
  integer, parameter :: k_p_mn55__he4_cr52   = 1094
  integer, parameter :: k_fe56__p_mn55   = 1095
  integer, parameter :: k_fe56__he4_cr52   = 1096
  integer, parameter :: k_p_fe56__he4_mn53   = 1097
  integer, parameter :: k_fe46__mn46__weak__wc12   = 1098
  integer, parameter :: k_fe46__p_mn45   = 1099
  integer, parameter :: k_fe46__p_cr45__weak__wc12   = 1100
  integer, parameter :: k_fe46__he4_cr42   = 1101
  integer, parameter :: k_p_fe46__co47   = 1102
  integer, parameter :: k_he4_fe46__ni50   = 1103
  integer, parameter :: k_he4_fe46__p_co49   = 1104
  integer, parameter :: k_co49__fe49__weak__bqa_pos_   = 1105
  integer, parameter :: k_co49__p_fe48   = 1106
  integer, parameter :: k_co49__he4_mn45   = 1107
  integer, parameter :: k_p_co49__ni50   = 1108
  integer, parameter :: k_p_co49__he4_fe46   = 1109
  integer, parameter :: k_he4_co49__p_ni52   = 1110
  integer, parameter :: k_ni52__co52__weak__wc12   = 1111
  integer, parameter :: k_ni52__p_co51   = 1112
  integer, parameter :: k_ni52__p_fe51__weak__wc12   = 1113
  integer, parameter :: k_ni52__he4_fe48   = 1114
  integer, parameter :: k_p_ni52__he4_co49   = 1115
  integer, parameter :: k_mn44__cr44__weak__bqa_pos_   = 1116
  integer, parameter :: k_mn44__p_cr43   = 1117
  integer, parameter :: k_mn44__he4_v40   = 1118
  integer, parameter :: k_p_mn44__fe45   = 1119
  integer, parameter :: k_he4_mn44__co48   = 1120
  integer, parameter :: k_he4_mn44__p_fe47   = 1121
  integer, parameter :: k_fe47__mn47__weak__wc12   = 1122
  integer, parameter :: k_fe47__p_mn46   = 1123
  integer, parameter :: k_fe47__p_cr46__weak__wc12   = 1124
  integer, parameter :: k_fe47__he4_cr43   = 1125
  integer, parameter :: k_p_fe47__co48   = 1126
  integer, parameter :: k_he4_fe47__ni51   = 1127
  integer, parameter :: k_p_fe47__he4_mn44   = 1128
  integer, parameter :: k_he4_fe47__p_co50   = 1129
  integer, parameter :: k_co50__fe50__weak__wc12   = 1130
  integer, parameter :: k_co50__p_fe49   = 1131
  integer, parameter :: k_co50__p_mn49__weak__wc12   = 1132
  integer, parameter :: k_co50__he4_mn46   = 1133
  integer, parameter :: k_p_co50__ni51   = 1134
  integer, parameter :: k_p_co50__he4_fe47   = 1135
  integer, parameter :: k_he4_co50__p_ni53   = 1136
  integer, parameter :: k_ni53__co53__weak__wc12   = 1137
  integer, parameter :: k_ni53__p_co52   = 1138
  integer, parameter :: k_ni53__p_fe52__weak__wc12   = 1139
  integer, parameter :: k_ni53__he4_fe49   = 1140
  integer, parameter :: k_p_ni53__he4_co50   = 1141
  integer, parameter :: k_co47__fe47__weak__bqa_pos_   = 1142
  integer, parameter :: k_co47__p_fe46   = 1143
  integer, parameter :: k_p_co47__ni48   = 1144
  integer, parameter :: k_he4_co47__p_ni50   = 1145
  integer, parameter :: k_ni50__co50__weak__wc12   = 1146
  integer, parameter :: k_ni50__p_co49   = 1147
  integer, parameter :: k_ni50__p_fe49__weak__wc12   = 1148
  integer, parameter :: k_ni50__he4_fe46   = 1149
  integer, parameter :: k_p_ni50__he4_co47   = 1150
  integer, parameter :: k_fe45__mn45__weak__wc17   = 1151
  integer, parameter :: k_fe45__p_mn44   = 1152
  integer, parameter :: k_fe45__p_cr44__weak__wc17   = 1153
  integer, parameter :: k_fe45__p_p_v43__weak__wc12   = 1154
  integer, parameter :: k_he4_fe45__ni49   = 1155
  integer, parameter :: k_he4_fe45__p_co48   = 1156
  integer, parameter :: k_fe45__p_p_p_ti42__weak__wc12   = 1157
  integer, parameter :: k_co48__fe48__weak__bqa_pos_   = 1158
  integer, parameter :: k_co48__p_fe47   = 1159
  integer, parameter :: k_co48__he4_mn44   = 1160
  integer, parameter :: k_p_co48__ni49   = 1161
  integer, parameter :: k_p_co48__he4_fe45   = 1162
  integer, parameter :: k_he4_co48__p_ni51   = 1163
  integer, parameter :: k_ni51__co51__weak__wc17   = 1164
  integer, parameter :: k_ni51__p_co50   = 1165
  integer, parameter :: k_ni51__p_fe50__weak__wc12   = 1166
  integer, parameter :: k_ni51__he4_fe47   = 1167
  integer, parameter :: k_p_ni51__he4_co48   = 1168
  integer, parameter :: k_ni48__co48__weak__wc17   = 1169
  integer, parameter :: k_ni48__p_co47   = 1170
  integer, parameter :: k_ni49__co49__weak__wc12   = 1171
  integer, parameter :: k_ni49__p_co48   = 1172
  integer, parameter :: k_ni49__p_fe48__weak__wc12   = 1173
  integer, parameter :: k_ni49__he4_fe45   = 1174

  ! reactvec indices
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4
  integer, parameter :: i_dqweak      = 5
  integer, parameter :: i_epart       = 6

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), zion(:), bion(:)
  real(rt), allocatable, save :: nion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, bion, nion, mion, wion
#endif

  !$acc declare create(aion, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 2722
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init()

    implicit none

    integer :: i

    ! Allocate ion info arrays
    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(bion(nspec))
    allocate(nion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    spec_names(jp)   = "hydrogen-1"
    spec_names(jd)   = "hydrogen-2"
    spec_names(jhe3)   = "helium-3"
    spec_names(jhe4)   = "helium-4"
    spec_names(jli7)   = "lithium-7"
    spec_names(jbe7)   = "beryllium-7"
    spec_names(jbe8)   = "beryllium-8"
    spec_names(jb8)   = "boron-8"
    spec_names(jc12)   = "carbon-12"
    spec_names(jc13)   = "carbon-13"
    spec_names(jn12)   = "nitrogen-12"
    spec_names(jn13)   = "nitrogen-13"
    spec_names(jn14)   = "nitrogen-14"
    spec_names(jn15)   = "nitrogen-15"
    spec_names(jo14)   = "oxygen-14"
    spec_names(jo15)   = "oxygen-15"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jo17)   = "oxygen-17"
    spec_names(jo18)   = "oxygen-18"
    spec_names(jf16)   = "fluorine-16"
    spec_names(jf17)   = "fluorine-17"
    spec_names(jf18)   = "fluorine-18"
    spec_names(jf19)   = "fluorine-19"
    spec_names(jf20)   = "fluorine-20"
    spec_names(jne17)   = "neon-17"
    spec_names(jne18)   = "neon-18"
    spec_names(jne19)   = "neon-19"
    spec_names(jne20)   = "neon-20"
    spec_names(jne21)   = "neon-21"
    spec_names(jne22)   = "neon-22"
    spec_names(jna19)   = "sodium-19"
    spec_names(jna20)   = "sodium-20"
    spec_names(jna21)   = "sodium-21"
    spec_names(jna22)   = "sodium-22"
    spec_names(jna23)   = "sodium-23"
    spec_names(jna24)   = "sodium-24"
    spec_names(jmg21)   = "magnesium-21"
    spec_names(jmg22)   = "magnesium-22"
    spec_names(jmg23)   = "magnesium-23"
    spec_names(jmg24)   = "magnesium-24"
    spec_names(jmg25)   = "magnesium-25"
    spec_names(jmg26)   = "magnesium-26"
    spec_names(jal22)   = "aluminum-22"
    spec_names(jal23)   = "aluminum-23"
    spec_names(jal24)   = "aluminum-24"
    spec_names(jal25)   = "aluminum-25"
    spec_names(jal26)   = "aluminum-26"
    spec_names(jal27)   = "aluminum-27"
    spec_names(jal28)   = "aluminum-28"
    spec_names(jsi24)   = "silicon-24"
    spec_names(jsi25)   = "silicon-25"
    spec_names(jsi26)   = "silicon-26"
    spec_names(jsi27)   = "silicon-27"
    spec_names(jsi28)   = "silicon-28"
    spec_names(jsi29)   = "silicon-29"
    spec_names(jsi30)   = "silicon-30"
    spec_names(jp26)   = "phosphorus-26"
    spec_names(jp27)   = "phosphorus-27"
    spec_names(jp28)   = "phosphorus-28"
    spec_names(jp29)   = "phosphorus-29"
    spec_names(jp30)   = "phosphorus-30"
    spec_names(jp31)   = "phosphorus-31"
    spec_names(jp32)   = "phosphorus-32"
    spec_names(jp33)   = "phosphorus-33"
    spec_names(js28)   = "sulfur-28"
    spec_names(js29)   = "sulfur-29"
    spec_names(js30)   = "sulfur-30"
    spec_names(js31)   = "sulfur-31"
    spec_names(js32)   = "sulfur-32"
    spec_names(js33)   = "sulfur-33"
    spec_names(js34)   = "sulfur-34"
    spec_names(js35)   = "sulfur-35"
    spec_names(jcl29)   = "chlorine-29"
    spec_names(jcl30)   = "chlorine-30"
    spec_names(jcl31)   = "chlorine-31"
    spec_names(jcl32)   = "chlorine-32"
    spec_names(jcl33)   = "chlorine-33"
    spec_names(jcl34)   = "chlorine-34"
    spec_names(jcl35)   = "chlorine-35"
    spec_names(jcl36)   = "chlorine-36"
    spec_names(jcl37)   = "chlorine-37"
    spec_names(jar31)   = "argon-31"
    spec_names(jar32)   = "argon-32"
    spec_names(jar33)   = "argon-33"
    spec_names(jar34)   = "argon-34"
    spec_names(jar35)   = "argon-35"
    spec_names(jar36)   = "argon-36"
    spec_names(jar37)   = "argon-37"
    spec_names(jar38)   = "argon-38"
    spec_names(jar39)   = "argon-39"
    spec_names(jk33)   = "potassium-33"
    spec_names(jk34)   = "potassium-34"
    spec_names(jk35)   = "potassium-35"
    spec_names(jk36)   = "potassium-36"
    spec_names(jk37)   = "potassium-37"
    spec_names(jk38)   = "potassium-38"
    spec_names(jk39)   = "potassium-39"
    spec_names(jk40)   = "potassium-40"
    spec_names(jk41)   = "potassium-41"
    spec_names(jca34)   = "calcium-34"
    spec_names(jca35)   = "calcium-35"
    spec_names(jca36)   = "calcium-36"
    spec_names(jca37)   = "calcium-37"
    spec_names(jca38)   = "calcium-38"
    spec_names(jca39)   = "calcium-39"
    spec_names(jca40)   = "calcium-40"
    spec_names(jca41)   = "calcium-41"
    spec_names(jca42)   = "calcium-42"
    spec_names(jca43)   = "calcium-43"
    spec_names(jca44)   = "calcium-44"
    spec_names(jsc36)   = "scandium-36"
    spec_names(jsc37)   = "scandium-37"
    spec_names(jsc38)   = "scandium-38"
    spec_names(jsc39)   = "scandium-39"
    spec_names(jsc40)   = "scandium-40"
    spec_names(jsc41)   = "scandium-41"
    spec_names(jsc42)   = "scandium-42"
    spec_names(jsc43)   = "scandium-43"
    spec_names(jsc44)   = "scandium-44"
    spec_names(jsc45)   = "scandium-45"
    spec_names(jsc46)   = "scandium-46"
    spec_names(jti38)   = "titanium-38"
    spec_names(jti39)   = "titanium-39"
    spec_names(jti40)   = "titanium-40"
    spec_names(jti41)   = "titanium-41"
    spec_names(jti42)   = "titanium-42"
    spec_names(jti43)   = "titanium-43"
    spec_names(jti44)   = "titanium-44"
    spec_names(jti45)   = "titanium-45"
    spec_names(jti46)   = "titanium-46"
    spec_names(jti47)   = "titanium-47"
    spec_names(jti48)   = "titanium-48"
    spec_names(jv40)   = "vanadium-40"
    spec_names(jv41)   = "vanadium-41"
    spec_names(jv42)   = "vanadium-42"
    spec_names(jv43)   = "vanadium-43"
    spec_names(jv44)   = "vanadium-44"
    spec_names(jv45)   = "vanadium-45"
    spec_names(jv46)   = "vanadium-46"
    spec_names(jv47)   = "vanadium-47"
    spec_names(jv48)   = "vanadium-48"
    spec_names(jv49)   = "vanadium-49"
    spec_names(jv50)   = "vanadium-50"
    spec_names(jcr42)   = "chromium-42"
    spec_names(jcr43)   = "chromium-43"
    spec_names(jcr44)   = "chromium-44"
    spec_names(jcr45)   = "chromium-45"
    spec_names(jcr46)   = "chromium-46"
    spec_names(jcr47)   = "chromium-47"
    spec_names(jcr48)   = "chromium-48"
    spec_names(jcr49)   = "chromium-49"
    spec_names(jcr50)   = "chromium-50"
    spec_names(jcr51)   = "chromium-51"
    spec_names(jcr52)   = "chromium-52"
    spec_names(jmn44)   = "manganese-44"
    spec_names(jmn45)   = "manganese-45"
    spec_names(jmn46)   = "manganese-46"
    spec_names(jmn47)   = "manganese-47"
    spec_names(jmn48)   = "manganese-48"
    spec_names(jmn49)   = "manganese-49"
    spec_names(jmn50)   = "manganese-50"
    spec_names(jmn51)   = "manganese-51"
    spec_names(jmn52)   = "manganese-52"
    spec_names(jmn53)   = "manganese-53"
    spec_names(jmn54)   = "manganese-54"
    spec_names(jmn55)   = "manganese-55"
    spec_names(jfe45)   = "iron-45"
    spec_names(jfe46)   = "iron-46"
    spec_names(jfe47)   = "iron-47"
    spec_names(jfe48)   = "iron-48"
    spec_names(jfe49)   = "iron-49"
    spec_names(jfe50)   = "iron-50"
    spec_names(jfe51)   = "iron-51"
    spec_names(jfe52)   = "iron-52"
    spec_names(jfe53)   = "iron-53"
    spec_names(jfe54)   = "iron-54"
    spec_names(jfe55)   = "iron-55"
    spec_names(jfe56)   = "iron-56"
    spec_names(jco47)   = "cobalt-47"
    spec_names(jco48)   = "cobalt-48"
    spec_names(jco49)   = "cobalt-49"
    spec_names(jco50)   = "cobalt-50"
    spec_names(jco51)   = "cobalt-51"
    spec_names(jco52)   = "cobalt-52"
    spec_names(jco53)   = "cobalt-53"
    spec_names(jco54)   = "cobalt-54"
    spec_names(jco55)   = "cobalt-55"
    spec_names(jco56)   = "cobalt-56"
    spec_names(jni48)   = "nickel-48"
    spec_names(jni49)   = "nickel-49"
    spec_names(jni50)   = "nickel-50"
    spec_names(jni51)   = "nickel-51"
    spec_names(jni52)   = "nickel-52"
    spec_names(jni53)   = "nickel-53"
    spec_names(jni54)   = "nickel-54"
    spec_names(jni55)   = "nickel-55"
    spec_names(jni56)   = "nickel-56"

    short_spec_names(jp)   = "h1"
    short_spec_names(jd)   = "h2"
    short_spec_names(jhe3)   = "he3"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jli7)   = "li7"
    short_spec_names(jbe7)   = "be7"
    short_spec_names(jbe8)   = "be8"
    short_spec_names(jb8)   = "b8"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jc13)   = "c13"
    short_spec_names(jn12)   = "n12"
    short_spec_names(jn13)   = "n13"
    short_spec_names(jn14)   = "n14"
    short_spec_names(jn15)   = "n15"
    short_spec_names(jo14)   = "o14"
    short_spec_names(jo15)   = "o15"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jo17)   = "o17"
    short_spec_names(jo18)   = "o18"
    short_spec_names(jf16)   = "f16"
    short_spec_names(jf17)   = "f17"
    short_spec_names(jf18)   = "f18"
    short_spec_names(jf19)   = "f19"
    short_spec_names(jf20)   = "f20"
    short_spec_names(jne17)   = "ne17"
    short_spec_names(jne18)   = "ne18"
    short_spec_names(jne19)   = "ne19"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jne21)   = "ne21"
    short_spec_names(jne22)   = "ne22"
    short_spec_names(jna19)   = "na19"
    short_spec_names(jna20)   = "na20"
    short_spec_names(jna21)   = "na21"
    short_spec_names(jna22)   = "na22"
    short_spec_names(jna23)   = "na23"
    short_spec_names(jna24)   = "na24"
    short_spec_names(jmg21)   = "mg21"
    short_spec_names(jmg22)   = "mg22"
    short_spec_names(jmg23)   = "mg23"
    short_spec_names(jmg24)   = "mg24"
    short_spec_names(jmg25)   = "mg25"
    short_spec_names(jmg26)   = "mg26"
    short_spec_names(jal22)   = "al22"
    short_spec_names(jal23)   = "al23"
    short_spec_names(jal24)   = "al24"
    short_spec_names(jal25)   = "al25"
    short_spec_names(jal26)   = "al26"
    short_spec_names(jal27)   = "al27"
    short_spec_names(jal28)   = "al28"
    short_spec_names(jsi24)   = "si24"
    short_spec_names(jsi25)   = "si25"
    short_spec_names(jsi26)   = "si26"
    short_spec_names(jsi27)   = "si27"
    short_spec_names(jsi28)   = "si28"
    short_spec_names(jsi29)   = "si29"
    short_spec_names(jsi30)   = "si30"
    short_spec_names(jp26)   = "p26"
    short_spec_names(jp27)   = "p27"
    short_spec_names(jp28)   = "p28"
    short_spec_names(jp29)   = "p29"
    short_spec_names(jp30)   = "p30"
    short_spec_names(jp31)   = "p31"
    short_spec_names(jp32)   = "p32"
    short_spec_names(jp33)   = "p33"
    short_spec_names(js28)   = "s28"
    short_spec_names(js29)   = "s29"
    short_spec_names(js30)   = "s30"
    short_spec_names(js31)   = "s31"
    short_spec_names(js32)   = "s32"
    short_spec_names(js33)   = "s33"
    short_spec_names(js34)   = "s34"
    short_spec_names(js35)   = "s35"
    short_spec_names(jcl29)   = "cl29"
    short_spec_names(jcl30)   = "cl30"
    short_spec_names(jcl31)   = "cl31"
    short_spec_names(jcl32)   = "cl32"
    short_spec_names(jcl33)   = "cl33"
    short_spec_names(jcl34)   = "cl34"
    short_spec_names(jcl35)   = "cl35"
    short_spec_names(jcl36)   = "cl36"
    short_spec_names(jcl37)   = "cl37"
    short_spec_names(jar31)   = "ar31"
    short_spec_names(jar32)   = "ar32"
    short_spec_names(jar33)   = "ar33"
    short_spec_names(jar34)   = "ar34"
    short_spec_names(jar35)   = "ar35"
    short_spec_names(jar36)   = "ar36"
    short_spec_names(jar37)   = "ar37"
    short_spec_names(jar38)   = "ar38"
    short_spec_names(jar39)   = "ar39"
    short_spec_names(jk33)   = "k33"
    short_spec_names(jk34)   = "k34"
    short_spec_names(jk35)   = "k35"
    short_spec_names(jk36)   = "k36"
    short_spec_names(jk37)   = "k37"
    short_spec_names(jk38)   = "k38"
    short_spec_names(jk39)   = "k39"
    short_spec_names(jk40)   = "k40"
    short_spec_names(jk41)   = "k41"
    short_spec_names(jca34)   = "ca34"
    short_spec_names(jca35)   = "ca35"
    short_spec_names(jca36)   = "ca36"
    short_spec_names(jca37)   = "ca37"
    short_spec_names(jca38)   = "ca38"
    short_spec_names(jca39)   = "ca39"
    short_spec_names(jca40)   = "ca40"
    short_spec_names(jca41)   = "ca41"
    short_spec_names(jca42)   = "ca42"
    short_spec_names(jca43)   = "ca43"
    short_spec_names(jca44)   = "ca44"
    short_spec_names(jsc36)   = "sc36"
    short_spec_names(jsc37)   = "sc37"
    short_spec_names(jsc38)   = "sc38"
    short_spec_names(jsc39)   = "sc39"
    short_spec_names(jsc40)   = "sc40"
    short_spec_names(jsc41)   = "sc41"
    short_spec_names(jsc42)   = "sc42"
    short_spec_names(jsc43)   = "sc43"
    short_spec_names(jsc44)   = "sc44"
    short_spec_names(jsc45)   = "sc45"
    short_spec_names(jsc46)   = "sc46"
    short_spec_names(jti38)   = "ti38"
    short_spec_names(jti39)   = "ti39"
    short_spec_names(jti40)   = "ti40"
    short_spec_names(jti41)   = "ti41"
    short_spec_names(jti42)   = "ti42"
    short_spec_names(jti43)   = "ti43"
    short_spec_names(jti44)   = "ti44"
    short_spec_names(jti45)   = "ti45"
    short_spec_names(jti46)   = "ti46"
    short_spec_names(jti47)   = "ti47"
    short_spec_names(jti48)   = "ti48"
    short_spec_names(jv40)   = "v40"
    short_spec_names(jv41)   = "v41"
    short_spec_names(jv42)   = "v42"
    short_spec_names(jv43)   = "v43"
    short_spec_names(jv44)   = "v44"
    short_spec_names(jv45)   = "v45"
    short_spec_names(jv46)   = "v46"
    short_spec_names(jv47)   = "v47"
    short_spec_names(jv48)   = "v48"
    short_spec_names(jv49)   = "v49"
    short_spec_names(jv50)   = "v50"
    short_spec_names(jcr42)   = "cr42"
    short_spec_names(jcr43)   = "cr43"
    short_spec_names(jcr44)   = "cr44"
    short_spec_names(jcr45)   = "cr45"
    short_spec_names(jcr46)   = "cr46"
    short_spec_names(jcr47)   = "cr47"
    short_spec_names(jcr48)   = "cr48"
    short_spec_names(jcr49)   = "cr49"
    short_spec_names(jcr50)   = "cr50"
    short_spec_names(jcr51)   = "cr51"
    short_spec_names(jcr52)   = "cr52"
    short_spec_names(jmn44)   = "mn44"
    short_spec_names(jmn45)   = "mn45"
    short_spec_names(jmn46)   = "mn46"
    short_spec_names(jmn47)   = "mn47"
    short_spec_names(jmn48)   = "mn48"
    short_spec_names(jmn49)   = "mn49"
    short_spec_names(jmn50)   = "mn50"
    short_spec_names(jmn51)   = "mn51"
    short_spec_names(jmn52)   = "mn52"
    short_spec_names(jmn53)   = "mn53"
    short_spec_names(jmn54)   = "mn54"
    short_spec_names(jmn55)   = "mn55"
    short_spec_names(jfe45)   = "fe45"
    short_spec_names(jfe46)   = "fe46"
    short_spec_names(jfe47)   = "fe47"
    short_spec_names(jfe48)   = "fe48"
    short_spec_names(jfe49)   = "fe49"
    short_spec_names(jfe50)   = "fe50"
    short_spec_names(jfe51)   = "fe51"
    short_spec_names(jfe52)   = "fe52"
    short_spec_names(jfe53)   = "fe53"
    short_spec_names(jfe54)   = "fe54"
    short_spec_names(jfe55)   = "fe55"
    short_spec_names(jfe56)   = "fe56"
    short_spec_names(jco47)   = "co47"
    short_spec_names(jco48)   = "co48"
    short_spec_names(jco49)   = "co49"
    short_spec_names(jco50)   = "co50"
    short_spec_names(jco51)   = "co51"
    short_spec_names(jco52)   = "co52"
    short_spec_names(jco53)   = "co53"
    short_spec_names(jco54)   = "co54"
    short_spec_names(jco55)   = "co55"
    short_spec_names(jco56)   = "co56"
    short_spec_names(jni48)   = "ni48"
    short_spec_names(jni49)   = "ni49"
    short_spec_names(jni50)   = "ni50"
    short_spec_names(jni51)   = "ni51"
    short_spec_names(jni52)   = "ni52"
    short_spec_names(jni53)   = "ni53"
    short_spec_names(jni54)   = "ni54"
    short_spec_names(jni55)   = "ni55"
    short_spec_names(jni56)   = "ni56"

    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jd)   = 1.11228300000000d+00
    ebind_per_nucleon(jhe3)   = 2.57268000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jli7)   = 5.60643900000000d+00
    ebind_per_nucleon(jbe7)   = 5.37154800000000d+00
    ebind_per_nucleon(jbe8)   = 7.06243500000000d+00
    ebind_per_nucleon(jb8)   = 4.71715500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jc13)   = 7.46984900000000d+00
    ebind_per_nucleon(jn12)   = 6.17010900000000d+00
    ebind_per_nucleon(jn13)   = 7.23886300000000d+00
    ebind_per_nucleon(jn14)   = 7.47561400000000d+00
    ebind_per_nucleon(jn15)   = 7.69946000000000d+00
    ebind_per_nucleon(jo14)   = 7.05227800000000d+00
    ebind_per_nucleon(jo15)   = 7.46369200000000d+00
    ebind_per_nucleon(jo16)   = 7.97620600000000d+00
    ebind_per_nucleon(jo17)   = 7.75072800000000d+00
    ebind_per_nucleon(jo18)   = 7.76709700000000d+00
    ebind_per_nucleon(jf16)   = 6.96373100000000d+00
    ebind_per_nucleon(jf17)   = 7.54232800000000d+00
    ebind_per_nucleon(jf18)   = 7.63163800000000d+00
    ebind_per_nucleon(jf19)   = 7.77901800000000d+00
    ebind_per_nucleon(jf20)   = 7.72013400000000d+00
    ebind_per_nucleon(jne17)   = 6.64049900000000d+00
    ebind_per_nucleon(jne18)   = 7.34125700000000d+00
    ebind_per_nucleon(jne19)   = 7.56734300000000d+00
    ebind_per_nucleon(jne20)   = 8.03224000000000d+00
    ebind_per_nucleon(jne21)   = 7.97171300000000d+00
    ebind_per_nucleon(jne22)   = 8.08046500000000d+00
    ebind_per_nucleon(jna19)   = 6.93788500000000d+00
    ebind_per_nucleon(jna20)   = 7.29849600000000d+00
    ebind_per_nucleon(jna21)   = 7.76554700000000d+00
    ebind_per_nucleon(jna22)   = 7.91566700000000d+00
    ebind_per_nucleon(jna23)   = 8.11149300000000d+00
    ebind_per_nucleon(jna24)   = 8.06348800000000d+00
    ebind_per_nucleon(jmg21)   = 7.10503100000000d+00
    ebind_per_nucleon(jmg22)   = 7.66276100000000d+00
    ebind_per_nucleon(jmg23)   = 7.90111500000000d+00
    ebind_per_nucleon(jmg24)   = 8.26070900000000d+00
    ebind_per_nucleon(jmg25)   = 8.22350200000000d+00
    ebind_per_nucleon(jmg26)   = 8.33387000000000d+00
    ebind_per_nucleon(jal22)   = 6.78200000000000d+00
    ebind_per_nucleon(jal23)   = 7.33572700000000d+00
    ebind_per_nucleon(jal24)   = 7.64958200000000d+00
    ebind_per_nucleon(jal25)   = 8.02113600000000d+00
    ebind_per_nucleon(jal26)   = 8.14976500000000d+00
    ebind_per_nucleon(jal27)   = 8.33155300000000d+00
    ebind_per_nucleon(jal28)   = 8.30989400000000d+00
    ebind_per_nucleon(jsi24)   = 7.16723200000000d+00
    ebind_per_nucleon(jsi25)   = 7.48011000000000d+00
    ebind_per_nucleon(jsi26)   = 7.92470800000000d+00
    ebind_per_nucleon(jsi27)   = 8.12434100000000d+00
    ebind_per_nucleon(jsi28)   = 8.44774400000000d+00
    ebind_per_nucleon(jsi29)   = 8.44863500000000d+00
    ebind_per_nucleon(jsi30)   = 8.52065400000000d+00
    ebind_per_nucleon(jp26)   = 7.19800000000000d+00
    ebind_per_nucleon(jp27)   = 7.66343800000000d+00
    ebind_per_nucleon(jp28)   = 7.90747900000000d+00
    ebind_per_nucleon(jp29)   = 8.25123600000000d+00
    ebind_per_nucleon(jp30)   = 8.35350600000000d+00
    ebind_per_nucleon(jp31)   = 8.48116700000000d+00
    ebind_per_nucleon(jp32)   = 8.46412000000000d+00
    ebind_per_nucleon(jp33)   = 8.51380600000000d+00
    ebind_per_nucleon(js28)   = 7.47879000000000d+00
    ebind_per_nucleon(js29)   = 7.74852000000000d+00
    ebind_per_nucleon(js30)   = 8.12270700000000d+00
    ebind_per_nucleon(js31)   = 8.28180000000000d+00
    ebind_per_nucleon(js32)   = 8.49312900000000d+00
    ebind_per_nucleon(js33)   = 8.49763000000000d+00
    ebind_per_nucleon(js34)   = 8.58349800000000d+00
    ebind_per_nucleon(js35)   = 8.53785000000000d+00
    ebind_per_nucleon(jcl29)   = 7.15883200000000d+00
    ebind_per_nucleon(jcl30)   = 7.48000000000000d+00
    ebind_per_nucleon(jcl31)   = 7.86920900000000d+00
    ebind_per_nucleon(jcl32)   = 8.07240400000000d+00
    ebind_per_nucleon(jcl33)   = 8.30475500000000d+00
    ebind_per_nucleon(jcl34)   = 8.39897000000000d+00
    ebind_per_nucleon(jcl35)   = 8.52027800000000d+00
    ebind_per_nucleon(jcl36)   = 8.52193100000000d+00
    ebind_per_nucleon(jcl37)   = 8.57028100000000d+00
    ebind_per_nucleon(jar31)   = 7.25200000000000d+00
    ebind_per_nucleon(jar32)   = 7.70000800000000d+00
    ebind_per_nucleon(jar33)   = 7.92895500000000d+00
    ebind_per_nucleon(jar34)   = 8.19767200000000d+00
    ebind_per_nucleon(jar35)   = 8.32746100000000d+00
    ebind_per_nucleon(jar36)   = 8.51990900000000d+00
    ebind_per_nucleon(jar37)   = 8.52713900000000d+00
    ebind_per_nucleon(jar38)   = 8.61428000000000d+00
    ebind_per_nucleon(jar39)   = 8.56259800000000d+00
    ebind_per_nucleon(jk33)   = 7.40700000000000d+00
    ebind_per_nucleon(jk34)   = 7.67000000000000d+00
    ebind_per_nucleon(jk35)   = 7.96584000000000d+00
    ebind_per_nucleon(jk36)   = 8.14221900000000d+00
    ebind_per_nucleon(jk37)   = 8.33984700000000d+00
    ebind_per_nucleon(jk38)   = 8.43805800000000d+00
    ebind_per_nucleon(jk39)   = 8.55702500000000d+00
    ebind_per_nucleon(jk40)   = 8.53809000000000d+00
    ebind_per_nucleon(jk41)   = 8.57607200000000d+00
    ebind_per_nucleon(jca34)   = 7.20400000000000d+00
    ebind_per_nucleon(jca35)   = 7.48700000000000d+00
    ebind_per_nucleon(jca36)   = 7.81587900000000d+00
    ebind_per_nucleon(jca37)   = 8.00345600000000d+00
    ebind_per_nucleon(jca38)   = 8.24004300000000d+00
    ebind_per_nucleon(jca39)   = 8.36967000000000d+00
    ebind_per_nucleon(jca40)   = 8.55130300000000d+00
    ebind_per_nucleon(jca41)   = 8.54670600000000d+00
    ebind_per_nucleon(jca42)   = 8.61656300000000d+00
    ebind_per_nucleon(jca43)   = 8.60066300000000d+00
    ebind_per_nucleon(jca44)   = 8.65817500000000d+00
    ebind_per_nucleon(jsc36)   = 7.18900000000000d+00
    ebind_per_nucleon(jsc37)   = 7.53200000000000d+00
    ebind_per_nucleon(jsc38)   = 7.75100000000000d+00
    ebind_per_nucleon(jsc39)   = 8.01345600000000d+00
    ebind_per_nucleon(jsc40)   = 8.17366900000000d+00
    ebind_per_nucleon(jsc41)   = 8.36919800000000d+00
    ebind_per_nucleon(jsc42)   = 8.44493300000000d+00
    ebind_per_nucleon(jsc43)   = 8.53082500000000d+00
    ebind_per_nucleon(jsc44)   = 8.55737900000000d+00
    ebind_per_nucleon(jsc45)   = 8.61893100000000d+00
    ebind_per_nucleon(jsc46)   = 8.62201200000000d+00
    ebind_per_nucleon(jti38)   = 7.33200000000000d+00
    ebind_per_nucleon(jti39)   = 7.57400000000000d+00
    ebind_per_nucleon(jti40)   = 7.86228600000000d+00
    ebind_per_nucleon(jti41)   = 8.03438800000000d+00
    ebind_per_nucleon(jti42)   = 8.25924700000000d+00
    ebind_per_nucleon(jti43)   = 8.35293200000000d+00
    ebind_per_nucleon(jti44)   = 8.53352000000000d+00
    ebind_per_nucleon(jti45)   = 8.55572200000000d+00
    ebind_per_nucleon(jti46)   = 8.65645100000000d+00
    ebind_per_nucleon(jti47)   = 8.66122700000000d+00
    ebind_per_nucleon(jti48)   = 8.72300600000000d+00
    ebind_per_nucleon(jv40)   = 7.31700000000000d+00
    ebind_per_nucleon(jv41)   = 7.62500000000000d+00
    ebind_per_nucleon(jv42)   = 7.82400000000000d+00
    ebind_per_nucleon(jv43)   = 8.06951200000000d+00
    ebind_per_nucleon(jv44)   = 8.21046300000000d+00
    ebind_per_nucleon(jv45)   = 8.38002900000000d+00
    ebind_per_nucleon(jv46)   = 8.48613000000000d+00
    ebind_per_nucleon(jv47)   = 8.58222500000000d+00
    ebind_per_nucleon(jv48)   = 8.62306100000000d+00
    ebind_per_nucleon(jv49)   = 8.68290800000000d+00
    ebind_per_nucleon(jv50)   = 8.69591800000000d+00
    ebind_per_nucleon(jcr42)   = 7.46400000000000d+00
    ebind_per_nucleon(jcr43)   = 7.68000000000000d+00
    ebind_per_nucleon(jcr44)   = 7.94800000000000d+00
    ebind_per_nucleon(jcr45)   = 8.08772800000000d+00
    ebind_per_nucleon(jcr46)   = 8.30382300000000d+00
    ebind_per_nucleon(jcr47)   = 8.40719500000000d+00
    ebind_per_nucleon(jcr48)   = 8.57226900000000d+00
    ebind_per_nucleon(jcr49)   = 8.61329100000000d+00
    ebind_per_nucleon(jcr50)   = 8.70103200000000d+00
    ebind_per_nucleon(jcr51)   = 8.71200500000000d+00
    ebind_per_nucleon(jcr52)   = 8.77598900000000d+00
    ebind_per_nucleon(jmn44)   = 7.46700000000000d+00
    ebind_per_nucleon(jmn45)   = 7.75300000000000d+00
    ebind_per_nucleon(jmn46)   = 7.91900000000000d+00
    ebind_per_nucleon(jmn47)   = 8.13531100000000d+00
    ebind_per_nucleon(jmn48)   = 8.27418500000000d+00
    ebind_per_nucleon(jmn49)   = 8.43992900000000d+00
    ebind_per_nucleon(jmn50)   = 8.53269600000000d+00
    ebind_per_nucleon(jmn51)   = 8.63377200000000d+00
    ebind_per_nucleon(jmn52)   = 8.67032900000000d+00
    ebind_per_nucleon(jmn53)   = 8.73417500000000d+00
    ebind_per_nucleon(jmn54)   = 8.73796500000000d+00
    ebind_per_nucleon(jmn55)   = 8.76502200000000d+00
    ebind_per_nucleon(jfe45)   = 7.31300000000000d+00
    ebind_per_nucleon(jfe46)   = 7.60900000000000d+00
    ebind_per_nucleon(jfe47)   = 7.78500000000000d+00
    ebind_per_nucleon(jfe48)   = 8.02300000000000d+00
    ebind_per_nucleon(jfe49)   = 8.16131100000000d+00
    ebind_per_nucleon(jfe50)   = 8.35402600000000d+00
    ebind_per_nucleon(jfe51)   = 8.46075900000000d+00
    ebind_per_nucleon(jfe52)   = 8.60957400000000d+00
    ebind_per_nucleon(jfe53)   = 8.64879900000000d+00
    ebind_per_nucleon(jfe54)   = 8.73638200000000d+00
    ebind_per_nucleon(jfe55)   = 8.74659500000000d+00
    ebind_per_nucleon(jfe56)   = 8.79035400000000d+00
    ebind_per_nucleon(jco47)   = 7.40100000000000d+00
    ebind_per_nucleon(jco48)   = 7.60000000000000d+00
    ebind_per_nucleon(jco49)   = 7.84200000000000d+00
    ebind_per_nucleon(jco50)   = 8.00100000000000d+00
    ebind_per_nucleon(jco51)   = 8.19325400000000d+00
    ebind_per_nucleon(jco52)   = 8.32588600000000d+00
    ebind_per_nucleon(jco53)   = 8.47765800000000d+00
    ebind_per_nucleon(jco54)   = 8.56921700000000d+00
    ebind_per_nucleon(jco55)   = 8.66961800000000d+00
    ebind_per_nucleon(jco56)   = 8.69483600000000d+00
    ebind_per_nucleon(jni48)   = 7.26500000000000d+00
    ebind_per_nucleon(jni49)   = 7.45700000000000d+00
    ebind_per_nucleon(jni50)   = 7.71600000000000d+00
    ebind_per_nucleon(jni51)   = 7.87500000000000d+00
    ebind_per_nucleon(jni52)   = 8.07900000000000d+00
    ebind_per_nucleon(jni53)   = 8.21707400000000d+00
    ebind_per_nucleon(jni54)   = 8.39303200000000d+00
    ebind_per_nucleon(jni55)   = 8.49732000000000d+00
    ebind_per_nucleon(jni56)   = 8.64277900000000d+00

    aion(jp)   = 1.00000000000000d+00
    aion(jd)   = 2.00000000000000d+00
    aion(jhe3)   = 3.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jli7)   = 7.00000000000000d+00
    aion(jbe7)   = 7.00000000000000d+00
    aion(jbe8)   = 8.00000000000000d+00
    aion(jb8)   = 8.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jc13)   = 1.30000000000000d+01
    aion(jn12)   = 1.20000000000000d+01
    aion(jn13)   = 1.30000000000000d+01
    aion(jn14)   = 1.40000000000000d+01
    aion(jn15)   = 1.50000000000000d+01
    aion(jo14)   = 1.40000000000000d+01
    aion(jo15)   = 1.50000000000000d+01
    aion(jo16)   = 1.60000000000000d+01
    aion(jo17)   = 1.70000000000000d+01
    aion(jo18)   = 1.80000000000000d+01
    aion(jf16)   = 1.60000000000000d+01
    aion(jf17)   = 1.70000000000000d+01
    aion(jf18)   = 1.80000000000000d+01
    aion(jf19)   = 1.90000000000000d+01
    aion(jf20)   = 2.00000000000000d+01
    aion(jne17)   = 1.70000000000000d+01
    aion(jne18)   = 1.80000000000000d+01
    aion(jne19)   = 1.90000000000000d+01
    aion(jne20)   = 2.00000000000000d+01
    aion(jne21)   = 2.10000000000000d+01
    aion(jne22)   = 2.20000000000000d+01
    aion(jna19)   = 1.90000000000000d+01
    aion(jna20)   = 2.00000000000000d+01
    aion(jna21)   = 2.10000000000000d+01
    aion(jna22)   = 2.20000000000000d+01
    aion(jna23)   = 2.30000000000000d+01
    aion(jna24)   = 2.40000000000000d+01
    aion(jmg21)   = 2.10000000000000d+01
    aion(jmg22)   = 2.20000000000000d+01
    aion(jmg23)   = 2.30000000000000d+01
    aion(jmg24)   = 2.40000000000000d+01
    aion(jmg25)   = 2.50000000000000d+01
    aion(jmg26)   = 2.60000000000000d+01
    aion(jal22)   = 2.20000000000000d+01
    aion(jal23)   = 2.30000000000000d+01
    aion(jal24)   = 2.40000000000000d+01
    aion(jal25)   = 2.50000000000000d+01
    aion(jal26)   = 2.60000000000000d+01
    aion(jal27)   = 2.70000000000000d+01
    aion(jal28)   = 2.80000000000000d+01
    aion(jsi24)   = 2.40000000000000d+01
    aion(jsi25)   = 2.50000000000000d+01
    aion(jsi26)   = 2.60000000000000d+01
    aion(jsi27)   = 2.70000000000000d+01
    aion(jsi28)   = 2.80000000000000d+01
    aion(jsi29)   = 2.90000000000000d+01
    aion(jsi30)   = 3.00000000000000d+01
    aion(jp26)   = 2.60000000000000d+01
    aion(jp27)   = 2.70000000000000d+01
    aion(jp28)   = 2.80000000000000d+01
    aion(jp29)   = 2.90000000000000d+01
    aion(jp30)   = 3.00000000000000d+01
    aion(jp31)   = 3.10000000000000d+01
    aion(jp32)   = 3.20000000000000d+01
    aion(jp33)   = 3.30000000000000d+01
    aion(js28)   = 2.80000000000000d+01
    aion(js29)   = 2.90000000000000d+01
    aion(js30)   = 3.00000000000000d+01
    aion(js31)   = 3.10000000000000d+01
    aion(js32)   = 3.20000000000000d+01
    aion(js33)   = 3.30000000000000d+01
    aion(js34)   = 3.40000000000000d+01
    aion(js35)   = 3.50000000000000d+01
    aion(jcl29)   = 2.90000000000000d+01
    aion(jcl30)   = 3.00000000000000d+01
    aion(jcl31)   = 3.10000000000000d+01
    aion(jcl32)   = 3.20000000000000d+01
    aion(jcl33)   = 3.30000000000000d+01
    aion(jcl34)   = 3.40000000000000d+01
    aion(jcl35)   = 3.50000000000000d+01
    aion(jcl36)   = 3.60000000000000d+01
    aion(jcl37)   = 3.70000000000000d+01
    aion(jar31)   = 3.10000000000000d+01
    aion(jar32)   = 3.20000000000000d+01
    aion(jar33)   = 3.30000000000000d+01
    aion(jar34)   = 3.40000000000000d+01
    aion(jar35)   = 3.50000000000000d+01
    aion(jar36)   = 3.60000000000000d+01
    aion(jar37)   = 3.70000000000000d+01
    aion(jar38)   = 3.80000000000000d+01
    aion(jar39)   = 3.90000000000000d+01
    aion(jk33)   = 3.30000000000000d+01
    aion(jk34)   = 3.40000000000000d+01
    aion(jk35)   = 3.50000000000000d+01
    aion(jk36)   = 3.60000000000000d+01
    aion(jk37)   = 3.70000000000000d+01
    aion(jk38)   = 3.80000000000000d+01
    aion(jk39)   = 3.90000000000000d+01
    aion(jk40)   = 4.00000000000000d+01
    aion(jk41)   = 4.10000000000000d+01
    aion(jca34)   = 3.40000000000000d+01
    aion(jca35)   = 3.50000000000000d+01
    aion(jca36)   = 3.60000000000000d+01
    aion(jca37)   = 3.70000000000000d+01
    aion(jca38)   = 3.80000000000000d+01
    aion(jca39)   = 3.90000000000000d+01
    aion(jca40)   = 4.00000000000000d+01
    aion(jca41)   = 4.10000000000000d+01
    aion(jca42)   = 4.20000000000000d+01
    aion(jca43)   = 4.30000000000000d+01
    aion(jca44)   = 4.40000000000000d+01
    aion(jsc36)   = 3.60000000000000d+01
    aion(jsc37)   = 3.70000000000000d+01
    aion(jsc38)   = 3.80000000000000d+01
    aion(jsc39)   = 3.90000000000000d+01
    aion(jsc40)   = 4.00000000000000d+01
    aion(jsc41)   = 4.10000000000000d+01
    aion(jsc42)   = 4.20000000000000d+01
    aion(jsc43)   = 4.30000000000000d+01
    aion(jsc44)   = 4.40000000000000d+01
    aion(jsc45)   = 4.50000000000000d+01
    aion(jsc46)   = 4.60000000000000d+01
    aion(jti38)   = 3.80000000000000d+01
    aion(jti39)   = 3.90000000000000d+01
    aion(jti40)   = 4.00000000000000d+01
    aion(jti41)   = 4.10000000000000d+01
    aion(jti42)   = 4.20000000000000d+01
    aion(jti43)   = 4.30000000000000d+01
    aion(jti44)   = 4.40000000000000d+01
    aion(jti45)   = 4.50000000000000d+01
    aion(jti46)   = 4.60000000000000d+01
    aion(jti47)   = 4.70000000000000d+01
    aion(jti48)   = 4.80000000000000d+01
    aion(jv40)   = 4.00000000000000d+01
    aion(jv41)   = 4.10000000000000d+01
    aion(jv42)   = 4.20000000000000d+01
    aion(jv43)   = 4.30000000000000d+01
    aion(jv44)   = 4.40000000000000d+01
    aion(jv45)   = 4.50000000000000d+01
    aion(jv46)   = 4.60000000000000d+01
    aion(jv47)   = 4.70000000000000d+01
    aion(jv48)   = 4.80000000000000d+01
    aion(jv49)   = 4.90000000000000d+01
    aion(jv50)   = 5.00000000000000d+01
    aion(jcr42)   = 4.20000000000000d+01
    aion(jcr43)   = 4.30000000000000d+01
    aion(jcr44)   = 4.40000000000000d+01
    aion(jcr45)   = 4.50000000000000d+01
    aion(jcr46)   = 4.60000000000000d+01
    aion(jcr47)   = 4.70000000000000d+01
    aion(jcr48)   = 4.80000000000000d+01
    aion(jcr49)   = 4.90000000000000d+01
    aion(jcr50)   = 5.00000000000000d+01
    aion(jcr51)   = 5.10000000000000d+01
    aion(jcr52)   = 5.20000000000000d+01
    aion(jmn44)   = 4.40000000000000d+01
    aion(jmn45)   = 4.50000000000000d+01
    aion(jmn46)   = 4.60000000000000d+01
    aion(jmn47)   = 4.70000000000000d+01
    aion(jmn48)   = 4.80000000000000d+01
    aion(jmn49)   = 4.90000000000000d+01
    aion(jmn50)   = 5.00000000000000d+01
    aion(jmn51)   = 5.10000000000000d+01
    aion(jmn52)   = 5.20000000000000d+01
    aion(jmn53)   = 5.30000000000000d+01
    aion(jmn54)   = 5.40000000000000d+01
    aion(jmn55)   = 5.50000000000000d+01
    aion(jfe45)   = 4.50000000000000d+01
    aion(jfe46)   = 4.60000000000000d+01
    aion(jfe47)   = 4.70000000000000d+01
    aion(jfe48)   = 4.80000000000000d+01
    aion(jfe49)   = 4.90000000000000d+01
    aion(jfe50)   = 5.00000000000000d+01
    aion(jfe51)   = 5.10000000000000d+01
    aion(jfe52)   = 5.20000000000000d+01
    aion(jfe53)   = 5.30000000000000d+01
    aion(jfe54)   = 5.40000000000000d+01
    aion(jfe55)   = 5.50000000000000d+01
    aion(jfe56)   = 5.60000000000000d+01
    aion(jco47)   = 4.70000000000000d+01
    aion(jco48)   = 4.80000000000000d+01
    aion(jco49)   = 4.90000000000000d+01
    aion(jco50)   = 5.00000000000000d+01
    aion(jco51)   = 5.10000000000000d+01
    aion(jco52)   = 5.20000000000000d+01
    aion(jco53)   = 5.30000000000000d+01
    aion(jco54)   = 5.40000000000000d+01
    aion(jco55)   = 5.50000000000000d+01
    aion(jco56)   = 5.60000000000000d+01
    aion(jni48)   = 4.80000000000000d+01
    aion(jni49)   = 4.90000000000000d+01
    aion(jni50)   = 5.00000000000000d+01
    aion(jni51)   = 5.10000000000000d+01
    aion(jni52)   = 5.20000000000000d+01
    aion(jni53)   = 5.30000000000000d+01
    aion(jni54)   = 5.40000000000000d+01
    aion(jni55)   = 5.50000000000000d+01
    aion(jni56)   = 5.60000000000000d+01

    zion(jp)   = 1.00000000000000d+00
    zion(jd)   = 1.00000000000000d+00
    zion(jhe3)   = 2.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jli7)   = 3.00000000000000d+00
    zion(jbe7)   = 4.00000000000000d+00
    zion(jbe8)   = 4.00000000000000d+00
    zion(jb8)   = 5.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jc13)   = 6.00000000000000d+00
    zion(jn12)   = 7.00000000000000d+00
    zion(jn13)   = 7.00000000000000d+00
    zion(jn14)   = 7.00000000000000d+00
    zion(jn15)   = 7.00000000000000d+00
    zion(jo14)   = 8.00000000000000d+00
    zion(jo15)   = 8.00000000000000d+00
    zion(jo16)   = 8.00000000000000d+00
    zion(jo17)   = 8.00000000000000d+00
    zion(jo18)   = 8.00000000000000d+00
    zion(jf16)   = 9.00000000000000d+00
    zion(jf17)   = 9.00000000000000d+00
    zion(jf18)   = 9.00000000000000d+00
    zion(jf19)   = 9.00000000000000d+00
    zion(jf20)   = 9.00000000000000d+00
    zion(jne17)   = 1.00000000000000d+01
    zion(jne18)   = 1.00000000000000d+01
    zion(jne19)   = 1.00000000000000d+01
    zion(jne20)   = 1.00000000000000d+01
    zion(jne21)   = 1.00000000000000d+01
    zion(jne22)   = 1.00000000000000d+01
    zion(jna19)   = 1.10000000000000d+01
    zion(jna20)   = 1.10000000000000d+01
    zion(jna21)   = 1.10000000000000d+01
    zion(jna22)   = 1.10000000000000d+01
    zion(jna23)   = 1.10000000000000d+01
    zion(jna24)   = 1.10000000000000d+01
    zion(jmg21)   = 1.20000000000000d+01
    zion(jmg22)   = 1.20000000000000d+01
    zion(jmg23)   = 1.20000000000000d+01
    zion(jmg24)   = 1.20000000000000d+01
    zion(jmg25)   = 1.20000000000000d+01
    zion(jmg26)   = 1.20000000000000d+01
    zion(jal22)   = 1.30000000000000d+01
    zion(jal23)   = 1.30000000000000d+01
    zion(jal24)   = 1.30000000000000d+01
    zion(jal25)   = 1.30000000000000d+01
    zion(jal26)   = 1.30000000000000d+01
    zion(jal27)   = 1.30000000000000d+01
    zion(jal28)   = 1.30000000000000d+01
    zion(jsi24)   = 1.40000000000000d+01
    zion(jsi25)   = 1.40000000000000d+01
    zion(jsi26)   = 1.40000000000000d+01
    zion(jsi27)   = 1.40000000000000d+01
    zion(jsi28)   = 1.40000000000000d+01
    zion(jsi29)   = 1.40000000000000d+01
    zion(jsi30)   = 1.40000000000000d+01
    zion(jp26)   = 1.50000000000000d+01
    zion(jp27)   = 1.50000000000000d+01
    zion(jp28)   = 1.50000000000000d+01
    zion(jp29)   = 1.50000000000000d+01
    zion(jp30)   = 1.50000000000000d+01
    zion(jp31)   = 1.50000000000000d+01
    zion(jp32)   = 1.50000000000000d+01
    zion(jp33)   = 1.50000000000000d+01
    zion(js28)   = 1.60000000000000d+01
    zion(js29)   = 1.60000000000000d+01
    zion(js30)   = 1.60000000000000d+01
    zion(js31)   = 1.60000000000000d+01
    zion(js32)   = 1.60000000000000d+01
    zion(js33)   = 1.60000000000000d+01
    zion(js34)   = 1.60000000000000d+01
    zion(js35)   = 1.60000000000000d+01
    zion(jcl29)   = 1.70000000000000d+01
    zion(jcl30)   = 1.70000000000000d+01
    zion(jcl31)   = 1.70000000000000d+01
    zion(jcl32)   = 1.70000000000000d+01
    zion(jcl33)   = 1.70000000000000d+01
    zion(jcl34)   = 1.70000000000000d+01
    zion(jcl35)   = 1.70000000000000d+01
    zion(jcl36)   = 1.70000000000000d+01
    zion(jcl37)   = 1.70000000000000d+01
    zion(jar31)   = 1.80000000000000d+01
    zion(jar32)   = 1.80000000000000d+01
    zion(jar33)   = 1.80000000000000d+01
    zion(jar34)   = 1.80000000000000d+01
    zion(jar35)   = 1.80000000000000d+01
    zion(jar36)   = 1.80000000000000d+01
    zion(jar37)   = 1.80000000000000d+01
    zion(jar38)   = 1.80000000000000d+01
    zion(jar39)   = 1.80000000000000d+01
    zion(jk33)   = 1.90000000000000d+01
    zion(jk34)   = 1.90000000000000d+01
    zion(jk35)   = 1.90000000000000d+01
    zion(jk36)   = 1.90000000000000d+01
    zion(jk37)   = 1.90000000000000d+01
    zion(jk38)   = 1.90000000000000d+01
    zion(jk39)   = 1.90000000000000d+01
    zion(jk40)   = 1.90000000000000d+01
    zion(jk41)   = 1.90000000000000d+01
    zion(jca34)   = 2.00000000000000d+01
    zion(jca35)   = 2.00000000000000d+01
    zion(jca36)   = 2.00000000000000d+01
    zion(jca37)   = 2.00000000000000d+01
    zion(jca38)   = 2.00000000000000d+01
    zion(jca39)   = 2.00000000000000d+01
    zion(jca40)   = 2.00000000000000d+01
    zion(jca41)   = 2.00000000000000d+01
    zion(jca42)   = 2.00000000000000d+01
    zion(jca43)   = 2.00000000000000d+01
    zion(jca44)   = 2.00000000000000d+01
    zion(jsc36)   = 2.10000000000000d+01
    zion(jsc37)   = 2.10000000000000d+01
    zion(jsc38)   = 2.10000000000000d+01
    zion(jsc39)   = 2.10000000000000d+01
    zion(jsc40)   = 2.10000000000000d+01
    zion(jsc41)   = 2.10000000000000d+01
    zion(jsc42)   = 2.10000000000000d+01
    zion(jsc43)   = 2.10000000000000d+01
    zion(jsc44)   = 2.10000000000000d+01
    zion(jsc45)   = 2.10000000000000d+01
    zion(jsc46)   = 2.10000000000000d+01
    zion(jti38)   = 2.20000000000000d+01
    zion(jti39)   = 2.20000000000000d+01
    zion(jti40)   = 2.20000000000000d+01
    zion(jti41)   = 2.20000000000000d+01
    zion(jti42)   = 2.20000000000000d+01
    zion(jti43)   = 2.20000000000000d+01
    zion(jti44)   = 2.20000000000000d+01
    zion(jti45)   = 2.20000000000000d+01
    zion(jti46)   = 2.20000000000000d+01
    zion(jti47)   = 2.20000000000000d+01
    zion(jti48)   = 2.20000000000000d+01
    zion(jv40)   = 2.30000000000000d+01
    zion(jv41)   = 2.30000000000000d+01
    zion(jv42)   = 2.30000000000000d+01
    zion(jv43)   = 2.30000000000000d+01
    zion(jv44)   = 2.30000000000000d+01
    zion(jv45)   = 2.30000000000000d+01
    zion(jv46)   = 2.30000000000000d+01
    zion(jv47)   = 2.30000000000000d+01
    zion(jv48)   = 2.30000000000000d+01
    zion(jv49)   = 2.30000000000000d+01
    zion(jv50)   = 2.30000000000000d+01
    zion(jcr42)   = 2.40000000000000d+01
    zion(jcr43)   = 2.40000000000000d+01
    zion(jcr44)   = 2.40000000000000d+01
    zion(jcr45)   = 2.40000000000000d+01
    zion(jcr46)   = 2.40000000000000d+01
    zion(jcr47)   = 2.40000000000000d+01
    zion(jcr48)   = 2.40000000000000d+01
    zion(jcr49)   = 2.40000000000000d+01
    zion(jcr50)   = 2.40000000000000d+01
    zion(jcr51)   = 2.40000000000000d+01
    zion(jcr52)   = 2.40000000000000d+01
    zion(jmn44)   = 2.50000000000000d+01
    zion(jmn45)   = 2.50000000000000d+01
    zion(jmn46)   = 2.50000000000000d+01
    zion(jmn47)   = 2.50000000000000d+01
    zion(jmn48)   = 2.50000000000000d+01
    zion(jmn49)   = 2.50000000000000d+01
    zion(jmn50)   = 2.50000000000000d+01
    zion(jmn51)   = 2.50000000000000d+01
    zion(jmn52)   = 2.50000000000000d+01
    zion(jmn53)   = 2.50000000000000d+01
    zion(jmn54)   = 2.50000000000000d+01
    zion(jmn55)   = 2.50000000000000d+01
    zion(jfe45)   = 2.60000000000000d+01
    zion(jfe46)   = 2.60000000000000d+01
    zion(jfe47)   = 2.60000000000000d+01
    zion(jfe48)   = 2.60000000000000d+01
    zion(jfe49)   = 2.60000000000000d+01
    zion(jfe50)   = 2.60000000000000d+01
    zion(jfe51)   = 2.60000000000000d+01
    zion(jfe52)   = 2.60000000000000d+01
    zion(jfe53)   = 2.60000000000000d+01
    zion(jfe54)   = 2.60000000000000d+01
    zion(jfe55)   = 2.60000000000000d+01
    zion(jfe56)   = 2.60000000000000d+01
    zion(jco47)   = 2.70000000000000d+01
    zion(jco48)   = 2.70000000000000d+01
    zion(jco49)   = 2.70000000000000d+01
    zion(jco50)   = 2.70000000000000d+01
    zion(jco51)   = 2.70000000000000d+01
    zion(jco52)   = 2.70000000000000d+01
    zion(jco53)   = 2.70000000000000d+01
    zion(jco54)   = 2.70000000000000d+01
    zion(jco55)   = 2.70000000000000d+01
    zion(jco56)   = 2.70000000000000d+01
    zion(jni48)   = 2.80000000000000d+01
    zion(jni49)   = 2.80000000000000d+01
    zion(jni50)   = 2.80000000000000d+01
    zion(jni51)   = 2.80000000000000d+01
    zion(jni52)   = 2.80000000000000d+01
    zion(jni53)   = 2.80000000000000d+01
    zion(jni54)   = 2.80000000000000d+01
    zion(jni55)   = 2.80000000000000d+01
    zion(jni56)   = 2.80000000000000d+01

    nion(jp)   = 0.00000000000000d+00
    nion(jd)   = 1.00000000000000d+00
    nion(jhe3)   = 1.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jli7)   = 4.00000000000000d+00
    nion(jbe7)   = 3.00000000000000d+00
    nion(jbe8)   = 4.00000000000000d+00
    nion(jb8)   = 3.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jc13)   = 7.00000000000000d+00
    nion(jn12)   = 5.00000000000000d+00
    nion(jn13)   = 6.00000000000000d+00
    nion(jn14)   = 7.00000000000000d+00
    nion(jn15)   = 8.00000000000000d+00
    nion(jo14)   = 6.00000000000000d+00
    nion(jo15)   = 7.00000000000000d+00
    nion(jo16)   = 8.00000000000000d+00
    nion(jo17)   = 9.00000000000000d+00
    nion(jo18)   = 1.00000000000000d+01
    nion(jf16)   = 7.00000000000000d+00
    nion(jf17)   = 8.00000000000000d+00
    nion(jf18)   = 9.00000000000000d+00
    nion(jf19)   = 1.00000000000000d+01
    nion(jf20)   = 1.10000000000000d+01
    nion(jne17)   = 7.00000000000000d+00
    nion(jne18)   = 8.00000000000000d+00
    nion(jne19)   = 9.00000000000000d+00
    nion(jne20)   = 1.00000000000000d+01
    nion(jne21)   = 1.10000000000000d+01
    nion(jne22)   = 1.20000000000000d+01
    nion(jna19)   = 8.00000000000000d+00
    nion(jna20)   = 9.00000000000000d+00
    nion(jna21)   = 1.00000000000000d+01
    nion(jna22)   = 1.10000000000000d+01
    nion(jna23)   = 1.20000000000000d+01
    nion(jna24)   = 1.30000000000000d+01
    nion(jmg21)   = 9.00000000000000d+00
    nion(jmg22)   = 1.00000000000000d+01
    nion(jmg23)   = 1.10000000000000d+01
    nion(jmg24)   = 1.20000000000000d+01
    nion(jmg25)   = 1.30000000000000d+01
    nion(jmg26)   = 1.40000000000000d+01
    nion(jal22)   = 9.00000000000000d+00
    nion(jal23)   = 1.00000000000000d+01
    nion(jal24)   = 1.10000000000000d+01
    nion(jal25)   = 1.20000000000000d+01
    nion(jal26)   = 1.30000000000000d+01
    nion(jal27)   = 1.40000000000000d+01
    nion(jal28)   = 1.50000000000000d+01
    nion(jsi24)   = 1.00000000000000d+01
    nion(jsi25)   = 1.10000000000000d+01
    nion(jsi26)   = 1.20000000000000d+01
    nion(jsi27)   = 1.30000000000000d+01
    nion(jsi28)   = 1.40000000000000d+01
    nion(jsi29)   = 1.50000000000000d+01
    nion(jsi30)   = 1.60000000000000d+01
    nion(jp26)   = 1.10000000000000d+01
    nion(jp27)   = 1.20000000000000d+01
    nion(jp28)   = 1.30000000000000d+01
    nion(jp29)   = 1.40000000000000d+01
    nion(jp30)   = 1.50000000000000d+01
    nion(jp31)   = 1.60000000000000d+01
    nion(jp32)   = 1.70000000000000d+01
    nion(jp33)   = 1.80000000000000d+01
    nion(js28)   = 1.20000000000000d+01
    nion(js29)   = 1.30000000000000d+01
    nion(js30)   = 1.40000000000000d+01
    nion(js31)   = 1.50000000000000d+01
    nion(js32)   = 1.60000000000000d+01
    nion(js33)   = 1.70000000000000d+01
    nion(js34)   = 1.80000000000000d+01
    nion(js35)   = 1.90000000000000d+01
    nion(jcl29)   = 1.20000000000000d+01
    nion(jcl30)   = 1.30000000000000d+01
    nion(jcl31)   = 1.40000000000000d+01
    nion(jcl32)   = 1.50000000000000d+01
    nion(jcl33)   = 1.60000000000000d+01
    nion(jcl34)   = 1.70000000000000d+01
    nion(jcl35)   = 1.80000000000000d+01
    nion(jcl36)   = 1.90000000000000d+01
    nion(jcl37)   = 2.00000000000000d+01
    nion(jar31)   = 1.30000000000000d+01
    nion(jar32)   = 1.40000000000000d+01
    nion(jar33)   = 1.50000000000000d+01
    nion(jar34)   = 1.60000000000000d+01
    nion(jar35)   = 1.70000000000000d+01
    nion(jar36)   = 1.80000000000000d+01
    nion(jar37)   = 1.90000000000000d+01
    nion(jar38)   = 2.00000000000000d+01
    nion(jar39)   = 2.10000000000000d+01
    nion(jk33)   = 1.40000000000000d+01
    nion(jk34)   = 1.50000000000000d+01
    nion(jk35)   = 1.60000000000000d+01
    nion(jk36)   = 1.70000000000000d+01
    nion(jk37)   = 1.80000000000000d+01
    nion(jk38)   = 1.90000000000000d+01
    nion(jk39)   = 2.00000000000000d+01
    nion(jk40)   = 2.10000000000000d+01
    nion(jk41)   = 2.20000000000000d+01
    nion(jca34)   = 1.40000000000000d+01
    nion(jca35)   = 1.50000000000000d+01
    nion(jca36)   = 1.60000000000000d+01
    nion(jca37)   = 1.70000000000000d+01
    nion(jca38)   = 1.80000000000000d+01
    nion(jca39)   = 1.90000000000000d+01
    nion(jca40)   = 2.00000000000000d+01
    nion(jca41)   = 2.10000000000000d+01
    nion(jca42)   = 2.20000000000000d+01
    nion(jca43)   = 2.30000000000000d+01
    nion(jca44)   = 2.40000000000000d+01
    nion(jsc36)   = 1.50000000000000d+01
    nion(jsc37)   = 1.60000000000000d+01
    nion(jsc38)   = 1.70000000000000d+01
    nion(jsc39)   = 1.80000000000000d+01
    nion(jsc40)   = 1.90000000000000d+01
    nion(jsc41)   = 2.00000000000000d+01
    nion(jsc42)   = 2.10000000000000d+01
    nion(jsc43)   = 2.20000000000000d+01
    nion(jsc44)   = 2.30000000000000d+01
    nion(jsc45)   = 2.40000000000000d+01
    nion(jsc46)   = 2.50000000000000d+01
    nion(jti38)   = 1.60000000000000d+01
    nion(jti39)   = 1.70000000000000d+01
    nion(jti40)   = 1.80000000000000d+01
    nion(jti41)   = 1.90000000000000d+01
    nion(jti42)   = 2.00000000000000d+01
    nion(jti43)   = 2.10000000000000d+01
    nion(jti44)   = 2.20000000000000d+01
    nion(jti45)   = 2.30000000000000d+01
    nion(jti46)   = 2.40000000000000d+01
    nion(jti47)   = 2.50000000000000d+01
    nion(jti48)   = 2.60000000000000d+01
    nion(jv40)   = 1.70000000000000d+01
    nion(jv41)   = 1.80000000000000d+01
    nion(jv42)   = 1.90000000000000d+01
    nion(jv43)   = 2.00000000000000d+01
    nion(jv44)   = 2.10000000000000d+01
    nion(jv45)   = 2.20000000000000d+01
    nion(jv46)   = 2.30000000000000d+01
    nion(jv47)   = 2.40000000000000d+01
    nion(jv48)   = 2.50000000000000d+01
    nion(jv49)   = 2.60000000000000d+01
    nion(jv50)   = 2.70000000000000d+01
    nion(jcr42)   = 1.80000000000000d+01
    nion(jcr43)   = 1.90000000000000d+01
    nion(jcr44)   = 2.00000000000000d+01
    nion(jcr45)   = 2.10000000000000d+01
    nion(jcr46)   = 2.20000000000000d+01
    nion(jcr47)   = 2.30000000000000d+01
    nion(jcr48)   = 2.40000000000000d+01
    nion(jcr49)   = 2.50000000000000d+01
    nion(jcr50)   = 2.60000000000000d+01
    nion(jcr51)   = 2.70000000000000d+01
    nion(jcr52)   = 2.80000000000000d+01
    nion(jmn44)   = 1.90000000000000d+01
    nion(jmn45)   = 2.00000000000000d+01
    nion(jmn46)   = 2.10000000000000d+01
    nion(jmn47)   = 2.20000000000000d+01
    nion(jmn48)   = 2.30000000000000d+01
    nion(jmn49)   = 2.40000000000000d+01
    nion(jmn50)   = 2.50000000000000d+01
    nion(jmn51)   = 2.60000000000000d+01
    nion(jmn52)   = 2.70000000000000d+01
    nion(jmn53)   = 2.80000000000000d+01
    nion(jmn54)   = 2.90000000000000d+01
    nion(jmn55)   = 3.00000000000000d+01
    nion(jfe45)   = 1.90000000000000d+01
    nion(jfe46)   = 2.00000000000000d+01
    nion(jfe47)   = 2.10000000000000d+01
    nion(jfe48)   = 2.20000000000000d+01
    nion(jfe49)   = 2.30000000000000d+01
    nion(jfe50)   = 2.40000000000000d+01
    nion(jfe51)   = 2.50000000000000d+01
    nion(jfe52)   = 2.60000000000000d+01
    nion(jfe53)   = 2.70000000000000d+01
    nion(jfe54)   = 2.80000000000000d+01
    nion(jfe55)   = 2.90000000000000d+01
    nion(jfe56)   = 3.00000000000000d+01
    nion(jco47)   = 2.00000000000000d+01
    nion(jco48)   = 2.10000000000000d+01
    nion(jco49)   = 2.20000000000000d+01
    nion(jco50)   = 2.30000000000000d+01
    nion(jco51)   = 2.40000000000000d+01
    nion(jco52)   = 2.50000000000000d+01
    nion(jco53)   = 2.60000000000000d+01
    nion(jco54)   = 2.70000000000000d+01
    nion(jco55)   = 2.80000000000000d+01
    nion(jco56)   = 2.90000000000000d+01
    nion(jni48)   = 2.00000000000000d+01
    nion(jni49)   = 2.10000000000000d+01
    nion(jni50)   = 2.20000000000000d+01
    nion(jni51)   = 2.30000000000000d+01
    nion(jni52)   = 2.40000000000000d+01
    nion(jni53)   = 2.50000000000000d+01
    nion(jni54)   = 2.60000000000000d+01
    nion(jni55)   = 2.70000000000000d+01
    nion(jni56)   = 2.80000000000000d+01

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    !wion(:) = aion(:)

    !$acc update device(aion, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec_evolve + 3)) ! neq + 1

    csr_jac_col_index = [ &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      170, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      180, &
      181, &
      182, &
      183, &
      184, &
      185, &
      186, &
      187, &
      188, &
      189, &
      190, &
      191, &
      192, &
      193, &
      194, &
      195, &
      196, &
      197, &
      198, &
      1, &
      2, &
      3, &
      4, &
      6, &
      198, &
      1, &
      2, &
      3, &
      4, &
      6, &
      198, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      8, &
      9, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      170, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      180, &
      181, &
      182, &
      183, &
      184, &
      185, &
      186, &
      187, &
      188, &
      190, &
      191, &
      192, &
      193, &
      194, &
      195, &
      196, &
      197, &
      198, &
      1, &
      4, &
      5, &
      6, &
      198, &
      1, &
      2, &
      3, &
      4, &
      6, &
      8, &
      198, &
      7, &
      8, &
      198, &
      1, &
      6, &
      8, &
      198, &
      1, &
      4, &
      9, &
      11, &
      12, &
      14, &
      17, &
      28, &
      198, &
      1, &
      10, &
      12, &
      13, &
      198, &
      1, &
      4, &
      11, &
      16, &
      198, &
      1, &
      4, &
      9, &
      12, &
      15, &
      17, &
      198, &
      1, &
      4, &
      10, &
      13, &
      15, &
      16, &
      18, &
      22, &
      198, &
      1, &
      4, &
      9, &
      14, &
      16, &
      17, &
      19, &
      23, &
      198, &
      1, &
      4, &
      12, &
      15, &
      21, &
      26, &
      198, &
      1, &
      4, &
      11, &
      13, &
      16, &
      20, &
      22, &
      27, &
      198, &
      1, &
      4, &
      9, &
      12, &
      14, &
      17, &
      21, &
      23, &
      25, &
      28, &
      32, &
      198, &
      1, &
      4, &
      13, &
      18, &
      21, &
      22, &
      24, &
      29, &
      198, &
      1, &
      4, &
      14, &
      19, &
      22, &
      23, &
      198, &
      1, &
      4, &
      16, &
      20, &
      27, &
      32, &
      198, &
      1, &
      4, &
      15, &
      17, &
      21, &
      25, &
      26, &
      28, &
      33, &
      198, &
      1, &
      4, &
      13, &
      16, &
      18, &
      22, &
      26, &
      27, &
      29, &
      34, &
      198, &
      1, &
      4, &
      14, &
      17, &
      19, &
      23, &
      27, &
      28, &
      30, &
      35, &
      198, &
      1, &
      4, &
      18, &
      24, &
      29, &
      198, &
      1, &
      4, &
      25, &
      32, &
      37, &
      198, &
      1, &
      4, &
      15, &
      21, &
      26, &
      31, &
      33, &
      38, &
      43, &
      198, &
      1, &
      4, &
      16, &
      20, &
      22, &
      27, &
      31, &
      32, &
      34, &
      39, &
      198, &
      1, &
      4, &
      9, &
      17, &
      21, &
      23, &
      24, &
      28, &
      32, &
      33, &
      35, &
      37, &
      40, &
      43, &
      45, &
      198, &
      1, &
      4, &
      18, &
      22, &
      24, &
      29, &
      33, &
      34, &
      36, &
      41, &
      198, &
      1, &
      4, &
      19, &
      23, &
      30, &
      34, &
      35, &
      42, &
      198, &
      1, &
      4, &
      26, &
      31, &
      38, &
      44, &
      198, &
      1, &
      4, &
      20, &
      25, &
      27, &
      32, &
      37, &
      39, &
      45, &
      198, &
      1, &
      4, &
      21, &
      26, &
      28, &
      33, &
      37, &
      38, &
      40, &
      43, &
      46, &
      198, &
      1, &
      4, &
      22, &
      27, &
      29, &
      34, &
      38, &
      39, &
      41, &
      44, &
      47, &
      198, &
      1, &
      4, &
      23, &
      28, &
      30, &
      35, &
      39, &
      40, &
      42, &
      45, &
      48, &
      198, &
      1, &
      4, &
      24, &
      29, &
      36, &
      41, &
      49, &
      198, &
      1, &
      4, &
      25, &
      32, &
      37, &
      43, &
      45, &
      51, &
      198, &
      1, &
      4, &
      26, &
      31, &
      33, &
      38, &
      43, &
      44, &
      46, &
      52, &
      198, &
      1, &
      4, &
      27, &
      32, &
      34, &
      39, &
      44, &
      45, &
      47, &
      50, &
      53, &
      198, &
      1, &
      4, &
      28, &
      33, &
      35, &
      40, &
      45, &
      46, &
      48, &
      51, &
      54, &
      59, &
      198, &
      1, &
      4, &
      29, &
      34, &
      36, &
      41, &
      46, &
      47, &
      49, &
      55, &
      198, &
      1, &
      4, &
      30, &
      35, &
      42, &
      47, &
      48, &
      56, &
      198, &
      1, &
      4, &
      37, &
      43, &
      51, &
      57, &
      198, &
      1, &
      4, &
      31, &
      38, &
      44, &
      50, &
      52, &
      58, &
      198, &
      1, &
      4, &
      32, &
      37, &
      39, &
      45, &
      50, &
      51, &
      53, &
      59, &
      198, &
      1, &
      4, &
      33, &
      38, &
      40, &
      46, &
      51, &
      52, &
      54, &
      60, &
      198, &
      1, &
      4, &
      34, &
      39, &
      41, &
      47, &
      52, &
      53, &
      55, &
      58, &
      61, &
      198, &
      1, &
      4, &
      35, &
      40, &
      42, &
      48, &
      53, &
      54, &
      56, &
      59, &
      62, &
      198, &
      1, &
      4, &
      36, &
      41, &
      49, &
      55, &
      63, &
      198, &
      1, &
      4, &
      44, &
      50, &
      58, &
      65, &
      198, &
      1, &
      4, &
      37, &
      43, &
      45, &
      51, &
      57, &
      59, &
      66, &
      198, &
      1, &
      4, &
      38, &
      44, &
      46, &
      52, &
      57, &
      58, &
      60, &
      67, &
      198, &
      1, &
      4, &
      39, &
      45, &
      47, &
      53, &
      58, &
      59, &
      61, &
      65, &
      68, &
      198, &
      1, &
      4, &
      40, &
      46, &
      48, &
      54, &
      59, &
      60, &
      62, &
      66, &
      69, &
      76, &
      198, &
      1, &
      4, &
      41, &
      47, &
      49, &
      55, &
      60, &
      61, &
      63, &
      70, &
      198, &
      1, &
      4, &
      42, &
      48, &
      56, &
      61, &
      62, &
      64, &
      71, &
      198, &
      1, &
      4, &
      43, &
      51, &
      57, &
      66, &
      74, &
      198, &
      1, &
      4, &
      44, &
      50, &
      52, &
      58, &
      65, &
      67, &
      75, &
      198, &
      1, &
      4, &
      45, &
      51, &
      53, &
      59, &
      65, &
      66, &
      68, &
      76, &
      198, &
      1, &
      4, &
      46, &
      52, &
      54, &
      60, &
      66, &
      67, &
      69, &
      77, &
      82, &
      198, &
      1, &
      4, &
      47, &
      53, &
      55, &
      61, &
      67, &
      68, &
      70, &
      75, &
      78, &
      198, &
      1, &
      4, &
      48, &
      54, &
      56, &
      62, &
      68, &
      69, &
      71, &
      76, &
      79, &
      198, &
      1, &
      4, &
      49, &
      55, &
      63, &
      70, &
      72, &
      80, &
      198, &
      1, &
      4, &
      56, &
      64, &
      71, &
      81, &
      198, &
      1, &
      4, &
      50, &
      58, &
      65, &
      73, &
      75, &
      83, &
      198, &
      1, &
      4, &
      51, &
      57, &
      59, &
      66, &
      73, &
      74, &
      76, &
      84, &
      198, &
      1, &
      4, &
      52, &
      58, &
      60, &
      67, &
      74, &
      75, &
      77, &
      82, &
      85, &
      198, &
      1, &
      4, &
      53, &
      59, &
      61, &
      68, &
      75, &
      76, &
      78, &
      83, &
      86, &
      198, &
      1, &
      4, &
      54, &
      60, &
      62, &
      69, &
      76, &
      77, &
      79, &
      84, &
      87, &
      94, &
      198, &
      1, &
      4, &
      55, &
      61, &
      63, &
      70, &
      77, &
      78, &
      80, &
      88, &
      198, &
      1, &
      4, &
      56, &
      62, &
      64, &
      71, &
      78, &
      79, &
      81, &
      89, &
      198, &
      1, &
      4, &
      63, &
      72, &
      80, &
      90, &
      198, &
      1, &
      4, &
      65, &
      73, &
      83, &
      91, &
      198, &
      1, &
      4, &
      57, &
      66, &
      74, &
      82, &
      84, &
      92, &
      198, &
      1, &
      4, &
      58, &
      65, &
      67, &
      75, &
      82, &
      83, &
      85, &
      93, &
      198, &
      1, &
      4, &
      59, &
      66, &
      68, &
      76, &
      83, &
      84, &
      86, &
      94, &
      198, &
      1, &
      4, &
      60, &
      67, &
      69, &
      77, &
      84, &
      85, &
      87, &
      95, &
      101, &
      198, &
      1, &
      4, &
      61, &
      68, &
      70, &
      78, &
      85, &
      86, &
      88, &
      93, &
      96, &
      198, &
      1, &
      4, &
      62, &
      69, &
      71, &
      79, &
      86, &
      87, &
      89, &
      94, &
      97, &
      198, &
      1, &
      4, &
      63, &
      70, &
      72, &
      80, &
      88, &
      90, &
      98, &
      198, &
      1, &
      4, &
      64, &
      71, &
      81, &
      88, &
      89, &
      99, &
      198, &
      1, &
      4, &
      74, &
      82, &
      92, &
      101, &
      198, &
      1, &
      4, &
      65, &
      73, &
      75, &
      83, &
      91, &
      93, &
      102, &
      198, &
      1, &
      4, &
      66, &
      74, &
      76, &
      84, &
      91, &
      92, &
      94, &
      103, &
      198, &
      1, &
      4, &
      67, &
      75, &
      77, &
      85, &
      92, &
      93, &
      95, &
      101, &
      104, &
      198, &
      1, &
      4, &
      68, &
      76, &
      78, &
      86, &
      93, &
      94, &
      96, &
      102, &
      105, &
      198, &
      1, &
      4, &
      69, &
      77, &
      79, &
      87, &
      94, &
      95, &
      97, &
      103, &
      106, &
      115, &
      198, &
      1, &
      4, &
      70, &
      78, &
      80, &
      88, &
      95, &
      96, &
      98, &
      107, &
      198, &
      1, &
      4, &
      71, &
      79, &
      81, &
      89, &
      96, &
      97, &
      99, &
      108, &
      198, &
      1, &
      4, &
      72, &
      80, &
      90, &
      98, &
      109, &
      198, &
      1, &
      4, &
      73, &
      83, &
      91, &
      100, &
      102, &
      112, &
      198, &
      1, &
      4, &
      74, &
      82, &
      84, &
      92, &
      100, &
      101, &
      103, &
      113, &
      198, &
      1, &
      4, &
      75, &
      83, &
      85, &
      93, &
      101, &
      102, &
      104, &
      114, &
      198, &
      1, &
      4, &
      76, &
      84, &
      86, &
      94, &
      102, &
      103, &
      105, &
      115, &
      198, &
      1, &
      4, &
      77, &
      85, &
      87, &
      95, &
      103, &
      104, &
      106, &
      116, &
      198, &
      1, &
      4, &
      78, &
      86, &
      88, &
      96, &
      104, &
      105, &
      107, &
      117, &
      198, &
      1, &
      4, &
      79, &
      87, &
      89, &
      97, &
      105, &
      106, &
      108, &
      115, &
      118, &
      198, &
      1, &
      4, &
      80, &
      88, &
      90, &
      98, &
      107, &
      109, &
      119, &
      198, &
      1, &
      4, &
      81, &
      89, &
      99, &
      107, &
      108, &
      110, &
      120, &
      198, &
      1, &
      4, &
      91, &
      100, &
      112, &
      122, &
      198, &
      1, &
      4, &
      82, &
      92, &
      101, &
      111, &
      113, &
      123, &
      198, &
      1, &
      4, &
      83, &
      91, &
      93, &
      102, &
      104, &
      111, &
      112, &
      114, &
      124, &
      198, &
      1, &
      4, &
      84, &
      92, &
      94, &
      103, &
      112, &
      113, &
      115, &
      125, &
      198, &
      1, &
      4, &
      85, &
      93, &
      95, &
      102, &
      104, &
      113, &
      114, &
      116, &
      123, &
      126, &
      198, &
      1, &
      4, &
      86, &
      94, &
      96, &
      105, &
      114, &
      115, &
      117, &
      124, &
      127, &
      198, &
      1, &
      4, &
      87, &
      95, &
      97, &
      106, &
      115, &
      116, &
      118, &
      125, &
      128, &
      145, &
      198, &
      1, &
      4, &
      88, &
      96, &
      98, &
      107, &
      116, &
      117, &
      119, &
      129, &
      198, &
      1, &
      4, &
      89, &
      97, &
      99, &
      108, &
      117, &
      118, &
      120, &
      130, &
      198, &
      1, &
      4, &
      90, &
      98, &
      109, &
      118, &
      119, &
      121, &
      131, &
      198, &
      1, &
      4, &
      99, &
      110, &
      119, &
      120, &
      132, &
      198, &
      1, &
      4, &
      101, &
      111, &
      123, &
      133, &
      198, &
      1, &
      4, &
      91, &
      100, &
      102, &
      112, &
      122, &
      124, &
      134, &
      198, &
      1, &
      4, &
      92, &
      101, &
      103, &
      113, &
      122, &
      123, &
      125, &
      135, &
      198, &
      1, &
      4, &
      93, &
      102, &
      104, &
      114, &
      123, &
      124, &
      126, &
      136, &
      198, &
      1, &
      4, &
      94, &
      103, &
      105, &
      115, &
      124, &
      125, &
      127, &
      137, &
      198, &
      1, &
      4, &
      95, &
      104, &
      106, &
      116, &
      125, &
      126, &
      128, &
      138, &
      145, &
      198, &
      1, &
      4, &
      96, &
      105, &
      107, &
      117, &
      126, &
      127, &
      129, &
      139, &
      198, &
      1, &
      4, &
      97, &
      106, &
      108, &
      118, &
      127, &
      128, &
      130, &
      140, &
      198, &
      1, &
      4, &
      98, &
      107, &
      109, &
      119, &
      128, &
      129, &
      131, &
      141, &
      198, &
      1, &
      4, &
      99, &
      108, &
      110, &
      120, &
      129, &
      130, &
      132, &
      142, &
      198, &
      1, &
      4, &
      109, &
      121, &
      131, &
      143, &
      198, &
      1, &
      4, &
      100, &
      112, &
      122, &
      134, &
      144, &
      198, &
      1, &
      4, &
      101, &
      111, &
      113, &
      123, &
      133, &
      135, &
      145, &
      198, &
      1, &
      4, &
      102, &
      112, &
      114, &
      124, &
      133, &
      134, &
      136, &
      146, &
      198, &
      1, &
      4, &
      103, &
      113, &
      115, &
      125, &
      134, &
      135, &
      137, &
      144, &
      147, &
      198, &
      1, &
      4, &
      104, &
      114, &
      116, &
      126, &
      135, &
      136, &
      138, &
      145, &
      148, &
      167, &
      198, &
      1, &
      4, &
      105, &
      115, &
      117, &
      127, &
      136, &
      137, &
      139, &
      146, &
      149, &
      198, &
      1, &
      4, &
      106, &
      116, &
      118, &
      128, &
      137, &
      138, &
      140, &
      147, &
      150, &
      198, &
      1, &
      4, &
      107, &
      117, &
      119, &
      129, &
      138, &
      139, &
      141, &
      151, &
      198, &
      1, &
      4, &
      108, &
      118, &
      120, &
      130, &
      139, &
      140, &
      142, &
      152, &
      198, &
      1, &
      4, &
      109, &
      119, &
      121, &
      131, &
      140, &
      141, &
      143, &
      153, &
      198, &
      1, &
      4, &
      110, &
      120, &
      132, &
      141, &
      142, &
      154, &
      198, &
      1, &
      4, &
      111, &
      123, &
      133, &
      145, &
      155, &
      198, &
      1, &
      4, &
      112, &
      122, &
      124, &
      134, &
      144, &
      146, &
      156, &
      198, &
      1, &
      4, &
      113, &
      123, &
      125, &
      135, &
      144, &
      145, &
      147, &
      157, &
      198, &
      1, &
      4, &
      114, &
      124, &
      126, &
      136, &
      145, &
      146, &
      148, &
      158, &
      167, &
      198, &
      1, &
      4, &
      115, &
      125, &
      127, &
      137, &
      146, &
      147, &
      149, &
      159, &
      198, &
      1, &
      4, &
      116, &
      126, &
      128, &
      138, &
      147, &
      148, &
      150, &
      157, &
      160, &
      198, &
      1, &
      4, &
      117, &
      127, &
      129, &
      139, &
      148, &
      149, &
      151, &
      161, &
      198, &
      1, &
      4, &
      118, &
      128, &
      130, &
      140, &
      149, &
      150, &
      152, &
      159, &
      162, &
      198, &
      1, &
      4, &
      119, &
      129, &
      131, &
      141, &
      150, &
      151, &
      153, &
      163, &
      198, &
      1, &
      4, &
      120, &
      130, &
      132, &
      142, &
      151, &
      152, &
      154, &
      164, &
      198, &
      1, &
      4, &
      121, &
      131, &
      143, &
      153, &
      165, &
      198, &
      1, &
      4, &
      122, &
      134, &
      144, &
      156, &
      168, &
      198, &
      1, &
      4, &
      123, &
      133, &
      135, &
      145, &
      155, &
      157, &
      169, &
      198, &
      1, &
      4, &
      124, &
      134, &
      136, &
      146, &
      155, &
      156, &
      158, &
      167, &
      170, &
      198, &
      1, &
      4, &
      125, &
      135, &
      137, &
      147, &
      156, &
      157, &
      159, &
      168, &
      171, &
      198, &
      1, &
      4, &
      126, &
      136, &
      138, &
      148, &
      157, &
      158, &
      160, &
      169, &
      172, &
      198, &
      1, &
      4, &
      127, &
      137, &
      139, &
      149, &
      158, &
      159, &
      161, &
      170, &
      173, &
      198, &
      1, &
      4, &
      128, &
      138, &
      140, &
      150, &
      159, &
      160, &
      162, &
      171, &
      174, &
      198, &
      1, &
      4, &
      129, &
      139, &
      141, &
      151, &
      160, &
      161, &
      163, &
      175, &
      198, &
      1, &
      4, &
      130, &
      140, &
      142, &
      152, &
      161, &
      162, &
      164, &
      176, &
      198, &
      1, &
      4, &
      131, &
      141, &
      143, &
      153, &
      162, &
      163, &
      165, &
      177, &
      198, &
      1, &
      4, &
      132, &
      142, &
      154, &
      163, &
      164, &
      166, &
      178, &
      198, &
      1, &
      4, &
      133, &
      145, &
      155, &
      167, &
      169, &
      180, &
      198, &
      1, &
      4, &
      134, &
      144, &
      146, &
      156, &
      167, &
      168, &
      170, &
      181, &
      198, &
      1, &
      4, &
      135, &
      145, &
      147, &
      157, &
      168, &
      169, &
      171, &
      182, &
      198, &
      1, &
      4, &
      136, &
      146, &
      148, &
      158, &
      169, &
      170, &
      172, &
      183, &
      198, &
      1, &
      4, &
      137, &
      147, &
      149, &
      159, &
      170, &
      171, &
      173, &
      184, &
      198, &
      1, &
      4, &
      138, &
      148, &
      150, &
      160, &
      171, &
      172, &
      174, &
      182, &
      185, &
      198, &
      1, &
      4, &
      139, &
      149, &
      151, &
      161, &
      172, &
      173, &
      175, &
      186, &
      198, &
      1, &
      4, &
      140, &
      150, &
      152, &
      162, &
      173, &
      174, &
      176, &
      187, &
      198, &
      1, &
      4, &
      141, &
      151, &
      153, &
      163, &
      174, &
      175, &
      177, &
      188, &
      198, &
      1, &
      4, &
      142, &
      152, &
      154, &
      164, &
      175, &
      176, &
      178, &
      198, &
      1, &
      4, &
      143, &
      153, &
      165, &
      177, &
      198, &
      1, &
      4, &
      154, &
      166, &
      177, &
      178, &
      198, &
      1, &
      4, &
      155, &
      167, &
      180, &
      190, &
      198, &
      1, &
      4, &
      144, &
      156, &
      168, &
      179, &
      181, &
      191, &
      198, &
      1, &
      4, &
      145, &
      155, &
      157, &
      169, &
      179, &
      180, &
      182, &
      192, &
      198, &
      1, &
      4, &
      146, &
      156, &
      158, &
      170, &
      180, &
      181, &
      183, &
      190, &
      193, &
      198, &
      1, &
      4, &
      147, &
      157, &
      159, &
      171, &
      181, &
      182, &
      184, &
      191, &
      194, &
      198, &
      1, &
      4, &
      148, &
      158, &
      160, &
      172, &
      182, &
      183, &
      185, &
      192, &
      195, &
      198, &
      1, &
      4, &
      149, &
      159, &
      161, &
      173, &
      183, &
      184, &
      186, &
      193, &
      196, &
      198, &
      1, &
      4, &
      150, &
      160, &
      162, &
      174, &
      184, &
      185, &
      187, &
      194, &
      197, &
      198, &
      1, &
      4, &
      151, &
      161, &
      163, &
      175, &
      185, &
      186, &
      188, &
      198, &
      1, &
      4, &
      152, &
      162, &
      164, &
      176, &
      186, &
      187, &
      198, &
      1, &
      4, &
      153, &
      163, &
      165, &
      177, &
      187, &
      188, &
      198, &
      1, &
      4, &
      154, &
      164, &
      166, &
      178, &
      188, &
      198, &
      1, &
      4, &
      168, &
      179, &
      189, &
      191, &
      198, &
      1, &
      4, &
      155, &
      167, &
      169, &
      180, &
      189, &
      190, &
      192, &
      198, &
      1, &
      4, &
      156, &
      168, &
      170, &
      181, &
      190, &
      191, &
      193, &
      198, &
      1, &
      4, &
      157, &
      169, &
      171, &
      182, &
      191, &
      192, &
      194, &
      198, &
      1, &
      4, &
      158, &
      170, &
      172, &
      183, &
      192, &
      193, &
      195, &
      198, &
      1, &
      4, &
      159, &
      171, &
      173, &
      184, &
      193, &
      194, &
      196, &
      198, &
      1, &
      4, &
      160, &
      172, &
      174, &
      185, &
      194, &
      195, &
      197, &
      198, &
      1, &
      4, &
      161, &
      173, &
      175, &
      186, &
      195, &
      196, &
      198, &
      1, &
      4, &
      162, &
      174, &
      176, &
      187, &
      196, &
      197, &
      198, &
      1, &
      4, &
      163, &
      175, &
      177, &
      188, &
      197, &
      198, &
      1, &
      179, &
      189, &
      198, &
      1, &
      4, &
      167, &
      180, &
      190, &
      198, &
      1, &
      4, &
      168, &
      179, &
      181, &
      191, &
      198, &
      1, &
      4, &
      169, &
      180, &
      182, &
      192, &
      198, &
      1, &
      4, &
      170, &
      181, &
      183, &
      193, &
      198, &
      1, &
      4, &
      171, &
      182, &
      184, &
      194, &
      198, &
      1, &
      4, &
      172, &
      183, &
      185, &
      195, &
      198, &
      1, &
      4, &
      173, &
      184, &
      186, &
      196, &
      198, &
      1, &
      4, &
      174, &
      185, &
      187, &
      197, &
      198, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      170, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      180, &
      181, &
      182, &
      183, &
      184, &
      185, &
      186, &
      187, &
      188, &
      189, &
      190, &
      191, &
      192, &
      193, &
      194, &
      195, &
      196, &
      197, &
      198, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      170, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      180, &
      181, &
      182, &
      183, &
      184, &
      185, &
      186, &
      187, &
      188, &
      189, &
      190, &
      191, &
      192, &
      193, &
      194, &
      195, &
      196, &
      197, &
      198, &
      199  ]

    csr_jac_row_count = [ &
      1, &
      198, &
      204, &
      210, &
      405, &
      410, &
      417, &
      420, &
      424, &
      433, &
      438, &
      443, &
      450, &
      459, &
      468, &
      475, &
      484, &
      496, &
      505, &
      512, &
      519, &
      529, &
      540, &
      551, &
      557, &
      563, &
      573, &
      584, &
      600, &
      611, &
      620, &
      627, &
      637, &
      649, &
      661, &
      673, &
      681, &
      690, &
      701, &
      713, &
      726, &
      737, &
      746, &
      753, &
      762, &
      773, &
      784, &
      796, &
      808, &
      816, &
      823, &
      833, &
      844, &
      856, &
      869, &
      880, &
      890, &
      898, &
      908, &
      919, &
      931, &
      943, &
      955, &
      964, &
      971, &
      980, &
      991, &
      1003, &
      1015, &
      1028, &
      1039, &
      1050, &
      1057, &
      1064, &
      1073, &
      1084, &
      1095, &
      1107, &
      1119, &
      1131, &
      1141, &
      1150, &
      1157, &
      1167, &
      1178, &
      1190, &
      1202, &
      1215, &
      1226, &
      1237, &
      1245, &
      1254, &
      1265, &
      1276, &
      1287, &
      1298, &
      1309, &
      1321, &
      1331, &
      1341, &
      1348, &
      1357, &
      1369, &
      1380, &
      1393, &
      1405, &
      1418, &
      1429, &
      1440, &
      1450, &
      1458, &
      1465, &
      1475, &
      1486, &
      1497, &
      1508, &
      1520, &
      1531, &
      1542, &
      1553, &
      1564, &
      1571, &
      1579, &
      1589, &
      1600, &
      1612, &
      1625, &
      1637, &
      1649, &
      1660, &
      1671, &
      1682, &
      1691, &
      1699, &
      1709, &
      1720, &
      1732, &
      1743, &
      1755, &
      1766, &
      1778, &
      1789, &
      1800, &
      1808, &
      1816, &
      1826, &
      1838, &
      1850, &
      1862, &
      1874, &
      1886, &
      1897, &
      1908, &
      1919, &
      1929, &
      1938, &
      1949, &
      1960, &
      1971, &
      1982, &
      1994, &
      2005, &
      2016, &
      2027, &
      2037, &
      2044, &
      2051, &
      2058, &
      2067, &
      2078, &
      2090, &
      2102, &
      2114, &
      2126, &
      2138, &
      2148, &
      2157, &
      2166, &
      2174, &
      2181, &
      2191, &
      2201, &
      2211, &
      2221, &
      2231, &
      2241, &
      2250, &
      2259, &
      2267, &
      2271, &
      2277, &
      2284, &
      2291, &
      2298, &
      2305, &
      2312, &
      2319, &
      2326, &
      2524, &
      2723  ]
#endif

  end subroutine actual_network_init


  subroutine actual_network_finalize()
    ! Deallocate storage arrays

    if (allocated(aion)) then
       deallocate(aion)
    endif

    if (allocated(zion)) then
       deallocate(zion)
    endif

    if (allocated(bion)) then
       deallocate(bion)
    endif

    if (allocated(nion)) then
       deallocate(nion)
    endif

    if (allocated(mion)) then
       deallocate(mion)
    endif

    if (allocated(wion)) then
       deallocate(wion)
    endif

#ifdef REACT_SPARSE_JACOBIAN
    if (allocated(csr_jac_col_index)) then
       deallocate(csr_jac_col_index)
    endif

    if (allocated(csr_jac_row_count)) then
       deallocate(csr_jac_row_count)
    endif
#endif

  end subroutine actual_network_finalize


  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate

    !$acc routine seq

    implicit none

    real(rt) :: dydt(nspec), enuc

    !$gpu

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
