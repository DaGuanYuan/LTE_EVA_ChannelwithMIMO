%LTE_7:��֯������������ѡ���о�
%%%%%%%%%%%��ʼ�������ź�����%%%%%%%%%%%%%%%%%%%%%%%%%%
num_bit = 128000;
Nsamp = 4;
M = 4;
k = log2(M);
Rb = 20*(10^6);    %bit rate 
Rs = 2*Rb / k;
Tb = 1/Rb;
Ts = 1/Rs;
rng('default');
rng(67);    %seed = 67
msg = randi([0,1],1,num_bit);  %generate messege
SNR_awgn  = 20000000;%awgn snr
snr = SNR_awgn;

%%%%%%%%%%%QPSK����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg = transpose(msg);
msg_syb = transpose(reshape(msg, size(msg, 2)*k, size(msg, 1)/k));     %IQ, get symbols
msg_de = bi2de(msg_syb);    %decimal symbols
grayencode = bitxor(0:M-1, floor((0:M-1)/2));   %generate gray list by order,[0,1,3,2]
msg_gry_enc = grayencode(msg_de+1);     %gray encode:0->0(The first),1->1(The second)��2->3(The third),3->2(The forth)��so by order,msg_enc+1
msg_qpsk = pskmod(msg_gry_enc,M, pi/4);    %QPSK modulation

%%%%%%%%%%%%%%%%%%%%%%%%%�����任%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_carrier = 1024;     %���ز�����
len_cp = 30;    %length of cp

%%%%%%%%%%��ӵ�Ƶ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pilotfrequency = 12;     %��Ƶ���
pilot_num = fix(length(msg_qpsk)/pilotfrequency);    %��Ƶ����
rest_num2 = pilotfrequency*(pilot_num + 1) - length(msg_qpsk);  %��Ƶ����
pilot_num =pilot_num + 1;
msg_qpsk(1, pilot_num * pilotfrequency) = 0;    %��0

msg_pilot = transpose(reshape(msg_qpsk, pilotfrequency,length(msg_qpsk)/pilotfrequency));   %trick, ��ӵ�Ƶ�������
msg_pilot_aft = zeros(length(msg_qpsk)/pilotfrequency,pilotfrequency+1);    %initialize
for i = 1:length(msg_qpsk)/pilotfrequency
   msg_pilot_aft(i,:) = [1, msg_pilot(i,:)] ;   %���뵼Ƶ
end
msg_pilot_aft_casd = reshape(transpose(msg_pilot_aft), 1, size(msg_pilot_aft,1)*size(msg_pilot_aft,2));     %ת����

%%%%%%%%%�����任%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_prl = fix(length(msg_pilot_aft_casd)/num_carrier);   %number of the groups which are divided in serial-to-parallel conversion
more_num3 = length(msg_pilot_aft_casd) - num_carrier*group_prl; %the same as interleaving
rest_num3 = num_carrier*(group_prl + 1) - length(msg_pilot_aft_casd);  %�ڶ��β�0����
rest_num = rest_num2 + rest_num3 ;   %�ܲ���ĸ���
group_prl = group_prl + 1;
msg_pilot_aft_casd(1, num_carrier * group_prl) = 1;   %��1��Ϊ�˰Ѳ��ϵĲ���Ҳ������뵼Ƶ����������ڲ�NaN

msg_prl = transpose(reshape(msg_pilot_aft_casd, num_carrier, group_prl));     %parallel signals
PLoc_1 = 1: (pilotfrequency+1) :length(msg_pilot_aft_casd);   %��Ƶλ��
DLoc_1 = setxor(1:length(msg_pilot_aft_casd), PLoc_1);      %����λ��

%%%%%%%%%%%%OFDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = sqrt(num_carrier);  %normalization factor
msg_idft = zeros(group_prl, num_carrier);
msg_ofdm = zeros(group_prl, num_carrier+len_cp);
msg_ray0 = zeros(num_carrier+len_cp, group_prl);
msg_ray = zeros(group_prl, num_carrier+len_cp);
msg_awgn = zeros(group_prl, num_carrier+len_cp);
msg_nocp = zeros(group_prl, num_carrier);
msg_nocp_ray = zeros(group_prl, num_carrier);
msg_nocp_est = zeros(group_prl, num_carrier);
msg_fft = zeros(group_prl, num_carrier);
MSG_PRL = zeros(group_prl, num_carrier);
MSG_NOCP_EST = zeros(group_prl, num_carrier);
H_RAY = zeros(group_prl, num_carrier);
% h_ray = zeros(group_prl, num_carrier+len_cp);
cp = zeros(group_prl, len_cp);  %initial
for j = 1 : group_prl
   msg_idft(j,:) = ifft(msg_prl(j,:))/N;    %idft
   cp(j,:) = msg_idft(j,num_carrier-len_cp+1:num_carrier);   %cp, 12bit
   msg_ofdm(j,:) = cat(2,cp(j,:),msg_idft(j,:));    %add cp,concatenate matrix
end
%%%%%%%%%%%%%%%Rayleigh_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_ofdm_re = reshape(transpose(msg_ofdm),1,size(msg_ofdm,1)*size(msg_ofdm,2));
rayleigh_chan = comm.RayleighChannel(...
   'SampleRate', Rs,...
   'PathDelays', [0,30,150,310,370,710,1090,1730,2510]*1e-9,...
   'AveragePathGains',[0,-1.5,-1.4,-3.6,-0.6,-9.1,-7.0,-12.0,-16.9],...
   'NormalizePathGains',true,...
   'MaximumDopplerShift',100,...
   'RandomStream', 'mt19937ar with seed',...
   'Seed', 110,...
   'PathGainsOutputPort',true);%rayleigh channel 
[msg_ray,PathGains1] = step(rayleigh_chan, transpose(msg_ofdm_re));     %signal after rayleigh channel   

%    rayleigh_chan = rayleighchan(Ts, 100, [0, 1.5e-5], [0, -2]);     %rayleigh channel
%    rayleigh_chan.StorePathGains = 1;
%    msg_ray(j,:) = filter(rayleigh_chan, msg_ofdm(j, :));    %signal after rayleigh channel
   %%%%%%%%%%%%%%%AWGN_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        msg_awgn(j,:) = awgn(msg_ofdm(j,:),SNR_awgn,'measured');     %not passing rayleigh
msg_awgn = awgn(msg_ray,SNR_awgn,'measured',19980611);      %after passing rayleigh
%        msg_awgn(j,:) = msg_ofdm(j,:);    %test, no noise
msg_awgn_re = transpose(reshape(msg_awgn, num_carrier+len_cp, length(msg_awgn)/(num_carrier+len_cp)));
   %%%%%%%%%%%%%%%%%%test_no awgn%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    msg_nocp(j,:) = msg_ray(j,13:num_carrier+len_cp);    %wipe off cp
%    msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft
for j = 1 : group_prl

   %%%%%%%%%%%%%%%DFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   msg_nocp(j,:) = msg_awgn_re(j,len_cp+1:num_carrier+len_cp);    %wipe off cp
   msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft


end

%%%%%%%%%%%%%%%%%%�����任%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_casd = reshape(transpose(msg_fft),1,size(msg_fft,1)*size(msg_fft,2));
%%%%%%%%%%%%%%%%%�ŵ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%LS channel estimate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = msg_casd(PLoc_1);     %notice that pilot is 1
msg_casd_len = length(msg_casd);
h_est = interp1(PLoc_1, h, 1:msg_casd_len);
h_est(msg_casd_len - rest_num + 1 :msg_casd_len) = [];   %ȥ������0
PLoc = 1: (pilotfrequency+1) :length(h_est); 
DLoc = setxor(1:length(h_est), PLoc);      %����λ��
h_ray = h_est(DLoc);


msg_casd(msg_casd_len - rest_num2 +1 : msg_casd_len) = [];     %ȥ������0
msg_casd = msg_casd(DLoc);
%%%%%%%%%%%%%%%%%%%�ŵ�����(f)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msg_casd_zf = msg_casd./h_ray;
%     msg_casd_zf_len = length(msg_casd_zf);
%     msg_casd_zf(msg_casd_zf_len - rest_num2 +1 : msg_casd_zf_len) = [];


%%%%%%%%%%%%%%%%%%QPSK���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_gry_demod = pskdemod(msg_casd_zf,M, pi/4);      %qpsk demodulate
[dummy, graydecode] = sort(grayencode); 
graydecode = graydecode - 1;
msg_demod = graydecode(msg_gry_demod+1)';
msg_demod_bi = de2bi(msg_demod,k)'; 
msg_demod_bi = msg_demod_bi(:);     %gray->normal


%%%%%%%%%%%%%%%%%%%����BER&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% msg(num_bit-34: num_bit) = [];      %for convolutional code
msg_wrong = abs(msg_demod_bi - msg);      %for turbo code
BER = sum(msg_wrong)/num_bit;


