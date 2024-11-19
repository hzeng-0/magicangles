function scattermult(A,dotsize)

[C,~,ic]=uniquetol(A, 0.001, 'ByRows', true);
count=accumarray(ic,1); 
set1=C(count==1,:); set2=C(count==2,:); set3=C(count>2,:);

scatter(set1(:,1),set1(:,2),dotsize,'blue','filled');
scatter(set2(:,1),set2(:,2),dotsize,'red','filled');
scatter(set3(:,1),set3(:,2),dotsize,'black','filled');

end