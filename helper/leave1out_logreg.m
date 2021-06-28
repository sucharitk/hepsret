function fit = leave1out_logreg(X,Y, method)

if ~exist('method', 'var'), method =1; end

[ntr,nfeat] = size(X);

indices = 1:ntr;

fit = NaN(ntr,1);
xx = NaN([ntr-1 nfeat ntr]);
yy= NaN(ntr-1,ntr);
x2 = NaN([1 nfeat ntr]);

if iscategorical(Y)
    fit = categorical(fit);
    yy = categorical(yy);
end

for nt = 1:ntr
    test = indices==nt;
    train = ~test;
    xx(:,:,nt) = X(train,:);
    yy(:,nt) = Y(train);
    x2(:,:,nt) = X(test,:);
    
end

parfor nt = 1:ntr
    fit(nt) = binary_classify(xx(:,:,nt), yy(:,nt), ...
        x2(:,:,nt), method);
end
% 
% switch method
%     case 1
%         parfor nt = 1:ntr
%                         
%             b = glmfit(xx(:,:,nt),yy(:,nt), ...
%                 'binomial','logit'); % Logistic regression
%             fit(nt) = round(glmval(b,x2(:,:,nt),'logit')');
%             
%         end
%     case 2
%         parfor nt = 1:ntr
%             
%             Mdl = fitclinear(xx(:,:,nt),yy(:,nt), 'Learner', 'logistic');
%             fit(nt) = predict(Mdl,x2(:,:,nt));
%             
%         end
%     case 3
%         parfor nt = 1:ntr
%             
%             Mdl = fitclinear(xx(:,:,nt),yy(:,nt),'Learner', 'svm');
%             fit(nt) = predict(Mdl,x2(:,:,nt));
%         end
%         
%     case 4
%         
%         % better do 1: 4 takes much longer
%         parfor nt = 1:ntr
%             Mdl = fitglm(xx(:,:,nt),yy(:,nt), 'link', 'logit',...
%                 'Distribution', 'binomial');
%             fit(nt) = round(predict(Mdl, x2(:,:,nt)));
%         end
%     case 5
%         
%         parfor nt = 1:ntr
%             
%             Mdl = fitcsvm(xx(:,:,nt), yy(:,nt));
%             fit(nt) = predict(Mdl,x2(:,:,nt));
%             
%         end
%         
%     case 6
%         
%         parfor nt = 1:ntr
%             
%             [B,FitInfo] = lassoglm(xx(:,:,nt),yy(:,nt),'binomial', ...
%                 'MaxIter',15, 'NumLambda',10);
%             [~,idxLambdaMinDeviance] = min(FitInfo.Deviance);
%             B0 = FitInfo.Intercept(idxLambdaMinDeviance);
%             coef = [B0; B(:,idxLambdaMinDeviance)];
%             
%             yhat = glmval(coef,x2(:,:,nt),'logit');
%             fit(nt) = round(yhat);
%         end
%         
%  
% end
end