

dataDir1='C:\data1'
dataDir2='C:\data2'
AnalysisDir='C:\analysis'

createP
p.show_output=0;
p.CandidateFindMethod=4;
p.GaussSwitch=3;

BFfactor=0.75;

%untreatFiles=dir('*.czi');
nucSum=[];
Allspots=[];
cellCounter=1;
for cond=1:2
    switch cond
        case 1
     dataDir=dataDir1;

            titleText='DMSO'
        case 2
    dataDir=  dataDir2;

            titleText='AB1'
    end
    cd(dataDir)

    untreatFiles=dir('*.czi');
%     untreatFiles=dir('*.lsm');
    for f=1:length(untreatFiles)
        clear imageData threeColour imageDataM imageDataTemp
        %% load data
        imageDataTemp(:,:,:)=imEx1(untreatFiles(f).name);
        %% reorder data 
         imageData=imageDataTemp;
             imageData(:,:,3)=imageDataTemp(:,:,3);
    imageData(:,:,2)=imageDataTemp(:,:,2);
      imageDataM(:,:,:)=mat2gray(  imageData(:,:,:));
%         imageData(:,:,3)=imageDataTemp(:,:,2); %DIC
%          imageData(:,:,1)=imageDataTemp(:,:,1); %Hoescht stain for nucleus
%          imageData(:,:,2)=imageDataTemp(:,:,1);
        %% make nice 3 colour image
        imageDataM(:,:,:)=mat2gray(  imageData(:,:,:));
        %      imageDataM(:,:,:)=imageData(:,:,:);
        imageDataM(:,:,3)=imageDataM(:,:,3)-0.2;
        threeColour(:,:,:)=cat(3,BFfactor*imageDataM(:,:,3),BFfactor*imageDataM(:,:,3),BFfactor*imageDataM(:,:,3));
        threeColour(:,:,2)= threeColour(:,:,2)+imageDataM(:,:,2);
        threeColour(:,:,3)= threeColour(:,:,3)+imageDataM(:,:,1);
        
        %% segment and crop
        mask=SimpleSegment3(imageData(:,:,1),300,0);
        figure;
        imshow(threeColour);
        hold on
        cellMask=segmentTrypanosomes3(mask, imageDataM, imageData);
        %loop over found nuclei
        for n=1:size(mask,3)
            cellCoord=[];
            nucProps=regionprops(mask(:,:,n),'centroid');
            % croppedImage=imcrop(threeColour,[nucProps.Centroid(1)-75,nucProps.Centroid(2)-75,150,150]);
            [colourCrop,rect2]=imcrop(threeColour,[nucProps.Centroid(1)-75,nucProps.Centroid(2)-75,149,149]);
            
            croppedImage1=imcrop(imageDataM(:,:,2),[nucProps.Centroid(1)-75,nucProps.Centroid(2)-75,149,149]);
            croppedImage2=imcrop(imageData(:,:,2),[nucProps.Centroid(1)-75,nucProps.Centroid(2)-75,149,149]);
            croppedMask=imcrop(mask(:,:,n),[nucProps.Centroid(1)-75,nucProps.Centroid(2)-75,149,149]);
            [SpotsCh1, ~, ~,~, ~, ~,~] = ADEMScode2_85(croppedImage2,p);
            [cellCoord(:,2), cellCoord(:,1)]=find(croppedMask);
            nucSum(1,cellCounter)=sum(croppedImage2(croppedMask==1));
                       nucSum(2,cellCounter)=cond-1;

            spotInd=ismember(round(SpotsCh1(:,1:2)),cellCoord,'rows');
            spots=SpotsCh1(spotInd,:);
            spots(:,13)=cond-1;
            spots(:,14)=cellCounter;
            Allspots=cat(1,Allspots,spots);
            [row,col]=find(bwperim(mask(:,:,n))==1);
                scatter(col,row,1,'m')
                            [row,col]=find(bwperim(cellMask(:,:,n))==1);
                scatter(col,row,0.5,'y')
                
                scatter(spots(:,1)+subplus(rect2(1)),spots(:,2)+subplus(rect2(2)),20,'xw')
                text(nucProps.Centroid(1),nucProps.Centroid(2),num2str(cellCounter),'color','y')
            cellCounter=cellCounter+1;

            
        end
        saveas(gcf,strcat(AnalysisDir, untreatFiles(f).name(1:end-3),'png'))
        
        
    end
   %save(strcat(analDir,'\fociTracking.mat'),'Allspots')
      xlswrite(AnalysisDir,Allspots(Allspots(:,13)==(cond-1),:),titleText)

   %xlswrite('E:\Google Drive\Data\Trypanosomes\Adam\New Analysis June 2018\KKT9.xlsx',Allspots(Allspots(:,13)==(cond-1),:),titleText)
end
subplot(3,2,1)
bar([mean(Allspots(Allspots(:,13)==0,5)),mean(Allspots(Allspots(:,13)==1,5))])
hold on
errorbar([1,2],[mean(Allspots(Allspots(:,13)==0,5)),mean(Allspots(Allspots(:,13)==1,5))],...
    [std(Allspots(Allspots(:,13)==0,5))/length(Allspots(Allspots(:,13)==0,5)).^0.5,std(Allspots(Allspots(:,13)==1,5))/length(Allspots(Allspots(:,13)==0,5)).^0.5],'.','LineWidth', 2,'color','k')
set(gca,'xtick',[1,2],'xticklabel',{'untreated','treated'})
ylabel('mean foci intensity')

if ttest2(Allspots(Allspots(:,13)==0,5),Allspots(Allspots(:,13)==1,5))==1
    text(1.95,1.2*(mean(Allspots(Allspots(:,13)==1,5))+std(Allspots(Allspots(:,13)==1,5))/length(Allspots(Allspots(:,13)==0,5)).^0.5),'*','FontSize',20)
end
subplot(3,2,2)

KDFplotColor(Allspots(Allspots(:,13)==0,5),0,[0,1,0]);
hold on
KDFplotColor(Allspots(Allspots(:,13)==1,5),0,[1,0,0]);
xlabel('foci intensity')
ylabel('probability')
legend('untreated','treated')
subplot(3,2,3)

bar([mean(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1)),mean(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))])
hold on
errorbar([1,2],[mean(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1)),mean(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))],...
    [std(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1))/length(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1)).^0.5,...
    std(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))/length(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1)).^0.5],'.','LineWidth', 2,'color','k')
set(gca,'xtick',[1,2],'xticklabel',{'untreated','treated'})
ylabel('mean number of foci')

if ttest2(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1),histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))==1
    text(1.95,1.2*(mean(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))+std(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))/length(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1)).^0.5),'*','FontSize',20)
end

ttest2(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1),histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1))
subplot(3,2,4)

histogram(histcounts(Allspots(Allspots(:,13)==1,14),'binwidth',1),'binwidth',1,'FaceColor',[1,0,0],'FaceAlpha',0.5,'Normalization','probability')
hold on
histogram(histcounts(Allspots(Allspots(:,13)==0,14),'binwidth',1),'binwidth',1,'FaceColor',[0,1,0],'FaceAlpha',0.5,'Normalization','probability')
xlabel('number of foci/cell')
ylabel('proportion of cells')
legend('treated','untreated')

subplot(3,2,5)
bar([mean(nucSum(1,nucSum(2,:)==0)),mean(nucSum(1,nucSum(2,:)==1))])
hold on
errorbar([1,2],[mean(nucSum(1,nucSum(2,:)==0)),mean(nucSum(1,nucSum(2,:)==1))],...
    [std(nucSum(1,nucSum(2,:)==0))/length(nucSum(1,nucSum(2,:)==0)).^0.5,...
    std(nucSum(1,nucSum(2,:)==1))/length(nucSum(1,nucSum(2,:)==1)).^0.5],'.','LineWidth', 2,'color','k')
set(gca,'xtick',[1,2],'xticklabel',{'untreated','treated'})
ylabel('mean nuclear intensity')

if ttest2(nucSum(1,nucSum(2,:)==0),nucSum(1,nucSum(2,:)==1))==1
    text(1.95,1.2*(mean(nucSum(1,nucSum(2,:)==1))+std(nucSum(1,nucSum(2,:)==1))/length(nucSum(1,nucSum(2,:)==1)).^0.5),'*','FontSize',20)
end

subplot(3,2,6)
KDFplotColor(nucSum(1,nucSum(2,:)==0),0,[0,1,0]);
hold on
KDFplotColor(nucSum(1,nucSum(2,:)==1),0,[1,0,0]);
xlabel('nuclear intensity')
ylabel('probability')
legend('untreated','treated')

% figure;
% KDFplot(Allspots(Allspots(:,13)==0,5));
% hold on
% KDFplot(Allspots(Allspots(:,13)==1,5));
% xlabel('Foci intensity')
% ylabel('Probability')
% legend('untreated','treated')
% 
% figure;
% subplot(1,2,1)
% 
% histogram(histcounts(Allspots(Allspots(:,13)==0,14),'BinWidth',1))
% xlabel('Foci/Cell')
% ylabel('No Cells')
% title('untreated')
% xlim([-1,4])
% subplot(1,2,2)
% 
% histogram(histcounts(Allspots(Allspots(:,13)==1,14),'BinWidth',1))
% xlabel('Foci/Cell')
% ylabel('No Cells')
% title('treated')
% xlim([-1,4])