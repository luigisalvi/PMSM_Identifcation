%N.B. Durante l'esecuzione del codice non chiudere le figure. 
Ceidx = 100:50:250;
widx = ["150", "300", "400", "500"];
vel_fit=[];
torque_fit=[];
datasets = strings(0,16);
%Loop creazione stringhe nome dataset
i = 1;
for j = 1:length(widx)
    for k = 1:length(Ceidx)
        datasets(i) = "V" + widx(j) + "_" + num2str(Ceidx(k));
        i = i + 1;
    end
end
%Loop main
for d = 1:length(datasets)
    disp("Now working on: " + datasets(d) + '(n°'+num2str(d)+')')
    fit=model_identification(datasets, d);
    vel_fit(:,d)= fit(1);
    torque_fit(:,d)= fit(2); 
end
close all 

%Plotting fit stats
fit_plot(length(datasets), vel_fit, torque_fit);

%% Functions 
function  fit= model_identification(datasets, dIndex)
%Creo una cartella per conservare le immagini
mkdir("FiguresNEW/" + datasets(dIndex));

%Caricare un file di condizioni operative
cd = load(datasets(dIndex) + ".mat"); %Current Dataset

%% Sistema da identificare
Ts=125;
idx = 5500:14500; %intervallo di regime
data = iddata([cd.w_mecc_rpm(idx) cd.Ce(idx)'], [cd.iaMis(idx) cd.ibMis(idx) cd.icMis(idx) cd.Cl(idx)], Ts); 
%advice(data)
data = detrend(data);
u = data.u;
y = data.y;
data.OutputName = {'Velocità [RPM]', 'Coppia [N/m]'};
%plot(data)
data.OutputName = ["w_m_e_c_c_R_P_M" "Ce" ];
data.InputName = ["iaMis" "ibMis" "icMis" "Cload"];

%% Splitting datasets
n = length(cd.Cm(idx)');
nt = round(n/2);
train = data(1:nt);
validation = data(nt+1:end);

%% Ordine ottimo 
ord_Max = 30;

 for k = 1:ord_Max
    m1 = arx(train, [k*eye(2) k*ones(2,4) ones(2,4)]); %modello ARX ad ordine k su dati di train
    

    %simulazione del modello 
    y_sim= sim(m1,u); 
    
    %errore di predizione
    ep = pe(m1,data); %prediction error
    err_pred = ep.OutputData;
    
    %errore di simulazione
    err_sim = y-y_sim;
    
    %indici di aderenza (varianza) dell'errore di predizione
    Jpred_train(k) = cov(err_pred(1:nt));
    Jpred_valid(k) = cov(err_pred(nt+1:n));
    
    %indici di aderenza (varianza) dell'errore di simulazione
    Jsim_train(k) = cov(err_sim(1:nt));
    Jsim_valid(k) = cov(err_sim(nt+1:n));
    
    
    %indici di ottimalità
    FPE(k) = fpe(m1); %final prediction error
    AIC(k) = aic(m1);%Akaike information criterion 
    MDL(k) = log(nt)*(k/nt) + log(Jpred_train(k)); %Minimum description length
    
end

%%Plotting
x=1:1:ord_Max;
 
% figure (1)
% subplot(2,1,1)
% plot(x,Jpred_train,  x, Jpred_valid);
% xlabel('n - numero di parametri'), ylabel('J predizione'),
% legend('J training set', 'J validation set')
% title('Andamento degli indici di aderenza in predizione e simulazione')
% subplot(2,1,2)
% plot(x, Jsim_train, x, Jsim_valid);
% xlabel('n - numero di parametri'), ylabel('J simulazione'),
% legend('J training set', 'J validation set')
% title_string1 = "Indici di aderenza | Dataset: "+ num2str(dIndex);
% title(title_string1);
% saveas(gcf, getImgPath(datasets(dIndex), "indici_aderenza.png"))
% 
% 
figure (2)
subplot(3,1,1)
plot(x,FPE);
xlabel('n - numero di parametri'), ylabel('FPE'),
title_string2 = "Indici di ottimalità | Dataset: "+ num2str(dIndex);
title(title_string2);
subplot(3,1,2)
plot(x, AIC);
xlabel('n - numero di parametri'), ylabel('AIC'),
subplot(3,1,3)
plot(x, MDL);
xlabel('n - numero di parametri'), ylabel('MDL'),
saveas(gcf, getImgPath(datasets(dIndex), "indice_ottimo.png"))

%% Modello lineare e compare
m = arx(train , [30*eye(2) 30*ones(2,4) ones(2,4)]);

[~,fit_cell,~]=compare(validation,m);
fit =fit_cell; %cell2mat(fit_cell);

figure(3)
compare(validation, m)
title_string3 = "Velocità-Coppia | Dataset: "+ num2str(dIndex);
title(title_string3);
saveas(gcf, getImgPath(datasets(dIndex), "model_comparison.png"))
saveas(gcf, getImgPath(datasets(dIndex), "model_comparison.fig"))

%AutoCorr-XCorr
figure(4)
resid(m,validation)
saveas(gcf, getImgPath(datasets(dIndex), "resid.png"))

end

function path = getImgPath(dataset, imgName)
    path = "FiguresNEW/" + dataset + "/" + imgName;
end
