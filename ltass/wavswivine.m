wavswivine = [2
4
5
6
8
9
12
20
26
28
30
31
33
37
40
41
43
47
48
51
52
53
55
58
61
64
65
68
69
70
71
74
75
76
79
80
83
89
90
92
94
95
96
97
98
99
100
101
107
110
112
113
114
115
116
118
119
121
122
124
126
127
128
129
130
132
133
134
135
136
142
144
147
148
149
151
155
156
158
160
161
163
168
171
172
174
176
179
181
186
187
188
191
192
193
194
196
201
202
203
204
205
206
208
211
214
217
219
222
223
224
226
227
228
234
235
237
238
240
243
245
246
249
250
251
254
256
257
259
262
266
268
275
276
279
285
287
288
290
292
293
295
298
299
300
302
305
306
307
308
316
321
324
325
332
333
334
338
340
341
348
349
355
359
360
361
362
363
367
371
374
375
376
379
383
386
387
388
390
392
400
401
405
407
408
409
411
412
413
414
415
416
418
419
420
421
423
424
425
427
431
432
433
434
436
437
444
446
448
450
452
453
454
456
457
462
463
470
471
472
473
474
475
479
483
484
485
486
487
488
490
491
492
495
498
499
500
501
502
503
506
507
508
509
511
512
513
514
516
517
522
523
528
530
531
536
538
540
542
549
554
556
557
558
559
561
562
564
567
568
569
570
571
575
576
578
579
582
583
584
585
586
588
589
591
592
596
597
601
606
607
610
611
612
613
615
617
618
619
620
621
623
624
625
627
631
632
636
640
641
642
643
645
646
647
648
649
650
651
652
654
655
659
660
661
665
668
669
670
672
675
678
684
687
689
690
691
696
697
700
701
703
704
707
708
710
714
715
721
725];

for i=1:length(wavswivine)
    wdzfilenames{i} = ['C:\temp\wivine\wdz' num2str(wavswivine(i)) '.wav'];
end