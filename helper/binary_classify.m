function Ytest = binary_classify(Xtrain, Ytrain, Xtest, method)
%
% different methods to perform binary classification depending on
% dimensionality of data or what type of method to use (logistic regression
% or svm) or whether to use penalised regression methods, etc.

%
%


switch method
    
    case 1
        % simple logistic regression
        
        %         parfor nt = 1:ntr
        
        b = glmfit(Xtrain, Ytrain, ...
            'binomial','logit'); % Logistic regression
        Ytest = round(glmval(b, Xtest, 'logit')');
        
        %         end
    case 2
        % logistic regression with matlab's method (a bit slower than the
        % first one but seems to not be overfitting with large number of
        % features) 
        
        %         parfor nt = 1:ntr
        
        Mdl = fitclinear(Xtrain, Ytrain,...
            'Learner', 'logistic');
        Ytest = predict(Mdl,Xtest);
        
        %         end
    case 3
        %         parfor nt = 1:ntr
        
        Mdl = fitclinear(Xtrain, Ytrain, ...
            'Learner', 'svm');
        Ytest = predict(Mdl,Xtest);
        %         end
        
    case 4
        
        % better do 1: 4 takes much longer
        %         parfor nt = 1:ntr
        Mdl = fitglm(Xtrain, Ytrain, ...
            'link', 'logit',...
            'Distribution', 'binomial');
        Ytest = round(predict(Mdl, Xtest));
        %         end
    case 5
        
        %         parfor nt = 1:ntr
        
        Mdl = fitcsvm(Xtrain, Ytrain);
        Ytest = predict(Mdl,Xtest);
        
        %         end
        
    case 6
        
        %         parfor nt = 1:ntr
        
        [B,FitInfo] = lassoglm(Xtrain, Ytrain, ...
            'binomial', ...
            'MaxIter',15, 'NumLambda',15,...
            'Alpha', 1);
        [~,idxLambdaMinDeviance] = min(FitInfo.Deviance);
        B0 = FitInfo.Intercept(idxLambdaMinDeviance);
        coef = [B0; B(:,idxLambdaMinDeviance)];
        
        Ytest = round(glmval(coef,Xtest,'logit'));
         
        %         end
        
        
end
end