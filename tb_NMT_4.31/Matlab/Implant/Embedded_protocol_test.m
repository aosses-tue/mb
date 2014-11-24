function result = Embedded_protocol_test

% Embedded_protocol_test: Test of Encode_embedded_protocol, Decode_embedded_protocol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tester(mfilename);

Gen_test_sequences;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f_mp1, idle_mp1] = Encode_embedded_protocol(q_mp1);

Open_log('mp1.log', 'new');
Log_sequence(f_mp1);
Log_sequence(idle_mp1);
expected_log = {
'   es   ms   as   ws   gs   ts'
'    2   24  102    2    8 5000'
'    4       104               '
'   es   ms   as   ws   gs   ts'
'    4   24    0    2    8 5000'
};
Tester(Compare_log_file(expected_log));
r_mp1 = Decode_embedded_protocol(f_mp1);
Tester(q_mp1, r_mp1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f_mp2, idle_mp2] = Encode_embedded_protocol(q_mp2);

Open_log('mp2.log', 'new');
Log_sequence(f_mp2);
Log_sequence(idle_mp2);
expected_log = {
'   es   ms   as   ws   gs   ts'
'   12   25  202    2    8 5000'
'   14       204               '
'   es   ms   as   ws   gs   ts'
'   14   25    0    2    8 5000'
};
Tester(Compare_log_file(expected_log));
r_mp2 = Decode_embedded_protocol(f_mp2);
Tester(q_mp2, r_mp2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f_mp12, idle_mp12] = Encode_embedded_protocol(q_mp12);

Open_log('mp12.log', 'new');
Log_sequence(f_mp12);
Log_sequence(idle_mp12);
expected_log = {
'   es   ms   as   ws   gs   ts'
'    3   30   33    2    8 5000'
'   13   30   43               '
'   24   25    0               '
'   es   ms   as   ws   gs   ts'
'   24   25    0    2    8 5000'
};
Tester(Compare_log_file(expected_log));
r_mp12 = Decode_embedded_protocol(f_mp12);
Tester(q_mp12, r_mp12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f1_mp12, idle1_mp12] = Encode_embedded_protocol(q1_mp12);

Open_log('q1_mp12.log', 'new');
Log_sequence(f1_mp12);
Log_sequence(idle1_mp12);
expected_log = {
'   es   ms   as   ws   gs   ts'
'   10   30  100    2    8  347'
'   es   ms   as   ws   gs   ts'
'   24   25    0    2    8  347'
};
Tester(Compare_log_file(expected_log));
r1_mp12 = Decode_embedded_protocol(f1_mp12);
Tester(q1_mp12, r1_mp12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fidle_mp12, idleidle_mp12] = Encode_embedded_protocol(qidle_mp12);

Open_log('qidle_mp12.log', 'new');
Log_sequence(fidle_mp12);
Log_sequence(idleidle_mp12);
expected_log = {
'   es   ms   as   ws   gs   ts'
'   24   25    0    2    8  347'
'   es   ms   as   ws   gs   ts'
'   24   25    0    2    8  347'
};
Tester(Compare_log_file(expected_log));
ridle_mp12 = Decode_embedded_protocol(fidle_mp12);
Tester(qidle_mp12, ridle_mp12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;	% Report result
