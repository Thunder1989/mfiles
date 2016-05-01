%plot kw_stat
acc =reshape(acc, 20,10)';
tpr =reshape(tpr, 20,10)';
fpr =reshape(fpr, 20,10)';
%%
figure
hold on
grid on
errorbar(1:20,mean(acc),std(acc),'r','LineWidth',2)
errorbar(1:20,mean(tpr),std(tpr),'b','LineWidth',2)
errorbar(1:20,mean(fpr),std(fpr),'k','LineWidth',2)
xlim([0 21])
title('ngram = 1,3')
legend('ACC','TPR','FPR','Location','SouthEast')