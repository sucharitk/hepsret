function tne = add_actichamp_S(nn)
% for actichamp the triggers have an S at the beginning, add the S and
% return the string

if nn<10
    tne = ['S  ' num2str(nn)];
elseif nn<100
    tne = ['S ' num2str(nn)];
else
    tne= ['S' num2str(nn)];
end

end