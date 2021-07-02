%BER of BPSK and QPSK.
EbNo = [2:2:8]; % 按实际测试调整，请注意折算
ber_BPSK =[0.2, 1e-2, 1e-4,1e-6];% 按实际测试调整，理论误码率：ber_BPSK =berawgn(EbNo,'psk',2,'nondiff' );
ber_QPSK =[0.2, 1e-2, 1e-4,1e-6];% 按实际测试调整，理论误码率：ber_QPSK = berawgn(EbNo,'psk',4,'nondiff' );

f1 = figure;
semilogy(EbNo,ber_BPSK,'r');
hold on;
semilogy(EbNo,ber_QPSK,'-.*b');
hold off;
xlabel('EB/N0(dB)');
ylabel('BER');
legend('BPSK','QPSK');
zoom on;

%%　理论误码率，供参考
% %BER of BPSK and QPSK
% EsNo = [2:2:10];
% ber_BPSK = berawgn(EsNo,'psk',2,'nondiff' );
% ber_QPSK = berawgn(EsNo-3,'psk',4,'nondiff' );
% 
% f2 = figure;
% semilogy(EsNo,ber_BPSK,'r');
% hold on;
% semilogy(EsNo,ber_QPSK,'-.*b');
% hold off;
% xlabel('ES/N0(dB)');
% ylabel('BER');
% legend('BPSK','QPSK');
% zoom on;
% grid on