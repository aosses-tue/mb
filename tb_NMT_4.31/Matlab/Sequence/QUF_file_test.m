function result = QUF_file_test

% QUF_file_test: Test of Save_sequence, Read_sequence for QUF files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Copyright: Cochlear Ltd
%      $Change: 86418 $
%    $Revision: #1 $
%    $DateTime: 2008/03/04 14:27:13 $
%      Authors: Brett Swanson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = Tester(mfilename);

Gen_test_sequences;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tester(Get_num_pulses(q1_mp12), 1);
Save_sequence(q1_mp12, 'q1_mp12');          % Save to file in working dir.
Tester(filecmp('q1_mp12.quf','s1.quf'));    % Compare to master file.
r1 = Read_sequence('q1_mp12.quf');          % Read file in working dir.
Tester(q1_mp12, r1);                        % Check that read matches original.

Tester(Get_num_pulses(qscan), 47104);
Save_sequence(qscan, 'qscan');              % Save to file in working dir.
Tester(filecmp('qscan.quf','scan.quf'));    % Compare to master file.
rscan = Read_sequence('qscan.quf');         % Read file in working dir.
Tester(qscan, rscan);                       % Check that read matches original.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if verbose;
    Disp_sequence(q1_mp12);
    disp(sprintf('Number of pulses = %d', Get_num_pulses(q1_mp12)));
    % Don't display entire qscan sequence, as it is too long.
    disp(sprintf('Number of pulses = %d', Get_num_pulses(qscan)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = Tester;    % Report result
