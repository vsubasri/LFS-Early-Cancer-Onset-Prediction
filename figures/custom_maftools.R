get_vcColors = function(alpha = 1, websafe = FALSE, named = TRUE){
  if(websafe){
    col = c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
            "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
            "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E",
            "#607D8B")
  }else{
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060', '#535c68')
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }
  
  if(named){
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','In_Frame_Del','Missense_Mutation','Silent','Nonsense_Mutation',
                           'RNA','Splice_Site','Intron','Frame_Shift_Ins','Frame_Shift_Del_Ins','Deletion','In_Frame_Ins',
                           'Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event', 'pathway')
  }
  
  col
}

get_titvCol = function(alpha = 1){
  col = c("#F44336", "#3F51B5", "#2196F3", "#4CAF50", "#FFC107", "#FF9800")
  #col = c('coral4', 'lightcyan4', 'cornflowerblue', 'lightsalmon1', 'forestgreen', 'deeppink3')
  col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  names(col) = c('C>T', 'C>G', 'C>A', 'T>A', 'T>C', 'T>G')
  col
}

validateMaf = function(maf, rdup = TRUE, isTCGA = isTCGA, chatty = TRUE){
  
  #necessary fields.
  required.fields = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
                      'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode')
  
  #Change column names to standard names; i.e, camel case
  for(i in 1:length(required.fields)){
    colId = suppressWarnings(grep(pattern = paste("^",required.fields[i],"$",sep=""), x = colnames(maf), ignore.case = TRUE))
    if(length(colId) > 0){
      colnames(maf)[colId] = required.fields[i]
    }
  }
  
  missing.fileds = required.fields[!required.fields %in% colnames(maf)] #check if any of them are missing
  
  if(length(missing.fileds) > 0){
    missing.fileds = paste(missing.fileds[1], sep = ',', collapse = ', ')
    stop(paste('missing required fields from MAF:', missing.fileds)) #stop if any of required.fields are missing
  }
  
  #convert "-" to "." in "Tumor_Sample_Barcode" to avoid complexity in naming
  #maf$Tumor_Sample_Barcode = gsub(pattern = '-', replacement = '.', x = as.character(maf$Tumor_Sample_Barcode))
  
  if(rdup){
    maf = maf[, variantId := paste(Chromosome, Start_Position, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2, sep = ':')]
    if(nrow(maf[duplicated(variantId)]) > 0){
      if(chatty){
        cat("--Removed",  nrow(maf[duplicated(variantId)]) ,"duplicated variants\n")
      }
      maf = maf[!duplicated(variantId)]
    }
    maf[,variantId := NULL]
  }
  
  if(nrow(maf[Hugo_Symbol %in% ""]) > 0){
    if(chatty){
      cat('--Found ', nrow(maf[Hugo_Symbol %in% ""]), ' variants with no Gene Symbols\n')
      #print(maf[Hugo_Symbol %in% "", required.fields, with = FALSE])
      cat("--Annotating them as 'UnknownGene' for convenience\n")
    }
    maf$Hugo_Symbol = ifelse(test = maf$Hugo_Symbol == "", yes = 'UnknownGene', no = maf$Hugo_Symbol)
  }
  
  if(nrow(maf[is.na(Hugo_Symbol)]) > 0){
    if(chatty){
      cat('--Found ', nrow(maf[is.na(Hugo_Symbol) > 0]), ' variants with no Gene Symbols\n')
      #print(maf[is.na(Hugo_Symbol), required.fields, with =FALSE])
      cat("--Annotating them as 'UnknownGene' for convenience\n")
    }
    maf$Hugo_Symbol = ifelse(test = is.na(maf$Hugo_Symbol), yes = 'UnknownGene', no = maf$Hugo_Symbol)
  }
  
  if(isTCGA){
    maf$Tumor_Sample_Barcode = substr(x = maf$Tumor_Sample_Barcode, start = 1, stop = 12)
  }
  
  #Variant Classification with Low/Modifier variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  silent = c("3'UTR", "5'UTR", "3'Flank", "Targeted_Region", "Silent", "Intron",
             "RNA", "IGR", "Splice_Region", "5'Flank", "lincRNA", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del")
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                   "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del","Frame_Shift_Del_Ins",
                   "In_Frame_Ins", "Missense_Mutation","Deletion")
  vt = c('SNP', 'DNP', 'TNP', 'ONP', 'INS', 'DEL','DEL_INS')
  
  maf.vcs = unique(as.character(maf[,Variant_Classification]))
  maf.vts = unique(as.character(maf[,Variant_Type]))
  
  if(length(maf.vcs[!maf.vcs %in% c(silent, vc.nonSilent)] > 0)){
    if(chatty){
      cat("--Non MAF specific values in Variant_Classification column:\n")
      for(temp in maf.vcs[!maf.vcs %in% c(silent, vc.nonSilent)]){
        cat(paste0("  ", temp, "\n"))
      }
    }
  }
  
  if(length(maf.vts[!maf.vts %in% vt] > 0)){
    if(chatty){
      cat("--Non MAF specific values in Variant_Type column:\n")
      for(temp in maf.vts[!maf.vts %in% vt]){
        cat(paste0("  ", temp, "\n"))
      }
    }
  }
  
  # Check type of variant position
  maf[,Chromosome := as.character(Chromosome)]
  maf[,Start_Position := as.numeric(as.character(Start_Position))]
  maf[,End_Position := as.numeric(as.character(End_Position))]
  #data.table::setkey(x = maf, Chromosome, Start_Position, End_Position)
  
  # Set Factors
  maf$Tumor_Sample_Barcode = as.factor(as.character(maf$Tumor_Sample_Barcode))
  
  return(maf)
}

dashboard = function(maf, color = NULL, rmOutlier = TRUE, log_conv = FALSE, titv.color = NULL, sfs = statFontSize, fontSize = fs, n = 10, donut = pie, rawcount = TRUE, stat = NULL, titleSize = NULL, barcodes = NULL, barcodeSize = NULL){
  
  if(is.null(color)){
    #hard coded color scheme if user doesnt provide any
    col = get_vcColors()
  }else{
    col = color
  }
  
  vcs = getSampleSummary(maf)
  vcs = vcs[,colnames(vcs)[!colnames(x = vcs) %in% c('total', 'Amp', 'Del', 'CNV_total')], with = FALSE]
  vcs = vcs[,c(1,order(colSums(x = vcs[,2:(ncol(vcs)), with =FALSE]), decreasing = TRUE)+1), with =FALSE] #order based on most event
  vcs.m = data.table::melt(data = vcs, id = 'Tumor_Sample_Barcode')
  colnames(vcs.m) = c('Tumor_Sample_Barcode', 'Variant_Classification', 'N')
  
  data.table::setDF(vcs)
  rownames(x = vcs) = vcs$Tumor_Sample_Barcode
  vcs = vcs[,-1, drop = FALSE]
  vcs = t(vcs)
  
  lo = matrix(data = 1:6, nrow = 2, byrow = TRUE)
  graphics::layout(mat = lo, heights = c(3.5, 3), widths = c(3, 2, 2))
  par(cex.axis = fontSize, font = 3, cex.main = titleSize[1], lwd = 1.2)
  
  #--------------------------- variant classification plot -----------------
  vc.plot.dat = rev(rowSums(vcs))
  if(log_conv){
    vc.plot.dat = log10(vc.plot.dat)
  }
  
  xt = pretty(c(0, vc.plot.dat))
  
  par(mar = c(3, 9, 3, 1))
  b = barplot(vc.plot.dat, axes = FALSE, horiz = TRUE, col = col[names(vc.plot.dat)], border = NA,
              xlim = c(0, max(xt)), names.arg = rep("", length(vc.plot.dat)))
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  axis(side = 2, at = b, labels = names(vc.plot.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.9, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Classification", adj = 0, cex.main = titleSize[1], font = 3)
  if(log_conv){
    axis(side = 2, at = 0, labels = "(log10)", lwd = 1.2, font = 3,
         las = 1, cex.axis = fontSize*0.9, hadj = 0.5, padj = 2, line = 0.75, tick = FALSE, outer = FALSE)
  }
  
  #--------------------------- variant type plot -----------------
  vt.plot.dat = maf@variant.type.summary
  vt.plot.dat = vt.plot.dat[,colnames(vt.plot.dat)[!colnames(x = vt.plot.dat) %in% c('total', 'CNV')], with = FALSE]
  vt.plot.dat = suppressWarnings(data.table::melt(vt.plot.dat[,c(2:(ncol(vt.plot.dat))), with = FALSE], id = NULL)[,sum(value), variable])
  colnames(vt.plot.dat)[2] = "sum"
  
  vt.cols = RColorBrewer::brewer.pal(n = 10, name = "Set3")
  if(log_conv){
    vt.plot.dat$sum = log10(vt.plot.dat$sum)
  }
  #xt = as.integer(seq(0, max(vt.plot.dat$sum), length.out = 4))
  xt = pretty(c(0, vt.plot.dat$sum))
  
  par(mar = c(3, 3, 3, 1))
  b = barplot(vt.plot.dat$sum, axes = FALSE, horiz = TRUE, col = vt.cols[1:length(vt.plot.dat$variable)],
              border = NA, xlim = c(0, max(xt)))
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  axis(side = 2, at = b, labels = vt.plot.dat$variable, lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "Variant Type", adj = 0, cex.main = titleSize[1], font = 3)
  if(log_conv){
    axis(side = 2, at = 0, labels = "(log10)", lwd = 1.2, font = 3,
         las = 1, cex.axis = fontSize*0.9, hadj = 0.5, padj = 2, line = 0.75, tick = FALSE, outer = FALSE)
  }
  
  #--------------------------- titv summary plot -----------------
  titv = titv(maf = maf, useSyn = TRUE, plot = FALSE)
  titv.counts = titv$raw.counts
  titv.sums = data.frame(value = colSums(titv.counts[,2:7]), stringsAsFactors = FALSE)
  titv.sums$class = rownames(titv.sums)
  if(!rawcount){
    titv.sums$raw_value = titv.sums$value
    titv.sums$value = titv.sums$value/sum(titv.sums$value)
    xt = seq(0, 1, 0.25)
  }else{
    xt = as.integer(seq(0, max(titv.sums$value, na.rm = TRUE), length.out = 4))
  }
  
  if(is.null(titv.color)){
    titv.color = get_titvCol()
  }
  
  par(mar = c(3, 3, 3, 1))
  b = barplot(titv.sums$value, axes = FALSE, horiz = TRUE, col = titv.color[rownames(titv.sums)],
              border = NA, xlim = c(0, xt[length(xt)]))
  axis(side = 2, at = b, labels = rownames(titv.sums), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = "SNV Class", adj = 0, cex.main = titleSize[1], font = 3)
  if(!rawcount){
    text(x = titv.sums$value+0.03, y = b, labels = titv.sums$raw_value,
         font = 4, col = "black", cex = fontSize, adj = 0)
  }
  abline(v = xt, h = 1:length(b)-0.25, lty = 2, lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  
  #--------------------------- variant per sample plot -----------------
  
  if(barcodes){
    par(mar = c(6, 2, 3, 1))
  } else{
    par(mar = c(3, 2, 3, 1))
  }
  
  if(log_conv){
    vcs = apply(vcs, 2, function(x){
      x_fract = x / sum(x)
      x_log_total = log10(sum(x))
      x_fract * x_log_total
    })
    #Replace NaN's and remove those samples (resulting from samples with no variants but with CNVs)
    vcs[is.nan(x = vcs)] = 0
    vcs = vcs[,-which(colSums(x = vcs) == 0)]
  }
  
  b = barplot(vcs, col = col[rownames(vcs)], border = NA, axes = FALSE, names.arg =  rep("", ncol(vcs)))
  axis(side = 2, at = as.integer(seq(0, max(colSums(vcs)), length.out = 4)), lwd = 1.2, font = 3, las = 2,
       line = -0.3, hadj = 0.6, cex.axis = fontSize)
  title(main = "Variants per sample", adj = 0, cex.main = titleSize[1], font = 3, line = 2)
  
  if(barcodes){
    mtext(text = colnames(vcs), side = 1, line = 0.2, outer = FALSE, las = 2, at = b, cex = barcodeSize)
  }
  
  if(!is.null(stat)){
    if(stat == 'mean'){
      med.line = round(maf@summary[ID %in% "total", Mean], 2)
      if(log_conv){
        med.line = round(log10(med.line), 2)
      }
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Mean: ', med.line, sep='')))
      title(main = paste0("Mean: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 3, line = 1)
      lines(x = c(1, b[length(b)]), y = c(med.line, med.line), col = "maroon", lwd = 1.2, lty = 2)
    }else if(stat == 'median'){
      med.line = round(maf@summary[ID %in% "total", Median], 2)
      if(log_conv){
        med.line = log10(med.line)
        title(main = paste0("Median: ", round(10^med.line)), adj = 0, cex.main = titleSize[1]*0.8, font = 3, line = 1)
      }else{
        title(main = paste0("Median: ", med.line), adj = 0, cex.main = titleSize[1]*0.8, font = 3, line = 1)
      }
      df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Median: ', med.line, sep='')))
      lines(x = c(1, b[length(b)]), y = c(med.line, med.line), col = "maroon", lwd = 1.2, lty = 2)
    }
  }
  
  if(log_conv){
    mtext(text = "(log10)", side = 2, lwd = 1.2, font = 3, cex = fontSize*0.7, line = 1)
  }
  
  #--------------------------- vc summary plot -----------------
  par(mar = c(3, 2, 3, 1))
  boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
  colnames(boxH)[ncol(boxH)] = 'boxStat'
  bcol = col[levels(vcs.m$Variant_Classification)]
  #vcs.m$N = log10(vcs.m$N)
  b = boxplot(N ~ Variant_Classification, data = vcs.m, xaxt="n", outline=FALSE, lty=1, lwd = 1.4, outwex=0,
              staplewex=0, axes = FALSE, border = bcol)
  
  axis(side = 2, at = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
       lwd = 1.2, font = 3, cex.axis = fontSize, las = 2)
  title(main = "Variant Classification \nsummary", adj = 0, cex.main = titleSize[1], font = 3, line = 1)
  abline(v = 1:length(bcol), h = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
         lty = 2,
         lwd = 0.3, col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
  
  #--------------------------- hugo-symbol plot -----------------
  gs = getGeneSummary(maf)
  nsamps = as.numeric(maf@summary[ID %in% "Samples", summary])
  gs.load = gs[,.(Hugo_Symbol, AlteredSamples)]
  gs.load[,AlteredSamples := round(AlteredSamples/nsamps, digits = 2) * 100]
  data.table::setDF(x = gs.load, rownames = gs.load$Hugo_Symbol)
  gs = gs[,colnames(gs)[!colnames(x = gs) %in% c('total', 'Amp', 'Del', 'CNV_total', 'MutatedSamples', 'AlteredSamples')], with = FALSE]
  
  if(nrow(gs) < n){
    gs.dat = gs
  }else{
    gs.dat = gs[1:n]
  }
  
  data.table::setDF(gs.dat)
  rownames(gs.dat) = gs.dat$Hugo_Symbol
  gs.dat = gs.dat[,-1, drop = FALSE]
  gs.dat = t(gs.dat)
  gs.dat = gs.dat[names(sort(rowSums(gs.dat), decreasing = TRUE)),, drop = FALSE]
  gs.dat = gs.dat[,names(sort(colSums(gs.dat))), drop = FALSE]
  
  xt = as.integer(seq(0, max(colSums(gs.dat))+2, length.out = 4))
  
  par(mar = c(3, 4, 3, 1))
  gs.load = gs.load[rev(colnames(gs.dat)),,]
  b = barplot(gs.dat, axes = FALSE, horiz = TRUE, col = col[rownames(gs.dat)], border = NA,
              xlim = c(0, max(xt)+(max(xt)*0.15)), names.arg = rep("", ncol(gs.dat)))
  axis(side = 2, at = b, labels = colnames(gs.dat), lwd = 1.2, cex.axis = fontSize,
       las = 2, line = 0.2, hadj = 0.8, font = 3, tick = FALSE)
  axis(side = 1, at = xt, lwd = 1.2, font = 3, las = 2, cex.axis = fontSize*0.9)
  title(main = paste0('Top ',  n, '\nmutated genes'), adj = 0, cex.main = titleSize[1], font = 3)
  text(x = colSums(gs.dat)+1, y = b, labels = rev(paste0(gs.load$AlteredSamples, "%")),
       font = 4, col = "black", cex = fontSize*0.9, adj = 0, xpd = TRUE)
  abline(h = b, v = xt,lty = 2, lwd = 0.3,
         col = grDevices::adjustcolor(col = "gray70", alpha.f = 0.6))
}

get_vcColors = function(alpha = 1, websafe = FALSE, named = TRUE){
  if(websafe){
    col = c("#F44336", "#E91E63", "#9C27B0", "#673AB7", "#3F51B5", "#2196F3",
            "#03A9F4", "#00BCD4", "#009688", "#4CAF50", "#8BC34A", "#CDDC39",
            "#FFEB3B", "#FFC107", "#FF9800", "#FF5722", "#795548", "#9E9E9E",
            "#607D8B")
  }else{
    col = c(RColorBrewer::brewer.pal(11, name = "Paired"), RColorBrewer::brewer.pal(11,name = "Spectral")[1:3],'black', 'violet', 'royalblue', '#7b7060', '#535c68')
    col = grDevices::adjustcolor(col = col, alpha.f = alpha)
  }
  
  if(named){
    names(col) = names = c('Nonstop_Mutation','Frame_Shift_Del','ITD','Missense_Mutation','Silent','Nonsense_Mutation',
                           'Splice_Site','Intron','Frame_Shift_Ins','In_Frame_Del','Frame_Shift_Del_Ins','In_Frame_Ins',
                           'Deletion','Translation_Start_Site',"Multi_Hit", 'Amp', 'Del', 'Complex_Event', 'pathway')
  }
  
  col
}

lollipopPlot_custom = function(maf, gene = NULL, AACol = NULL, labelPos = NULL, labPosSize = 0.9, showMutationRate = TRUE,
                               showDomainLabel = TRUE, cBioPortal = FALSE, refSeqID = NULL, proteinID = NULL, roundedRect = TRUE,
                               repel = FALSE, collapsePosLabel = TRUE, showLegend = TRUE, legendTxtSize = 0.8, labPosAngle = 0, domainLabelSize = 0.8, axisTextSize = c(1, 1),
                               printCount = FALSE, colors = NULL, domainAlpha = 1, domainBorderCol = "black", bgBorderCol = "black", labelOnlyUniqueDoamins = TRUE, defaultYaxis = FALSE, titleSize = c(1.2, 1), pointSize = 1.5){
  
  if(is.null(gene)){
    stop('Please provide a gene name.')
  }
  
  geneID = gene
  #Protein domain source.
  gff = system.file('extdata', 'protein_domains.RDs', package = 'maftools')
  gff = readRDS(file = gff)
  data.table::setDT(x = gff)
  
  mut = subsetMaf(maf = maf, includeSyn = FALSE, genes = gene, query = "Variant_Type != 'CNV'", mafObj = FALSE)
  
  if(is.null(AACol)){
    pchange = c('HGVSp_Short', 'Protein_Change', 'AAChange')
    if(length(pchange[pchange %in% colnames(mut)]) > 0){
      pchange = suppressWarnings(pchange[pchange %in% colnames(mut)][1])
      message(paste0("Assuming protein change information are stored under column ", pchange,". Use argument AACol to override if necessary."))
      colnames(mut)[which(colnames(mut) == pchange)] = 'AAChange_'
    }else{
      message('Available fields:')
      print(colnames(mut))
      stop('AAChange field not found in MAF. Use argument AACol to manually specifiy field name containing protein changes.')
    }
  }else{
    if(length(which(colnames(mut) == AACol)) == 0){
      message('Available fields:')
      print(colnames(mut))
      stop(paste0("Column ", AACol, " not found."))
    }else{
      colnames(mut)[which(colnames(mut) == AACol)] = 'AAChange_'
    }
  }
  
  prot.dat = mut[Hugo_Symbol %in% geneID, .(Variant_Type, Variant_Classification, AAChange_)]
  if(nrow(prot.dat) == 0){
    stop(paste(geneID, 'does not seem to have any mutations!', sep=' '))
  }
  
  prot = gff[HGNC %in% geneID]
  
  if(nrow(prot) == 0){
    stop(paste('Structure for protein', geneID, 'not found.', sep=' '))
  }
  
  if(!is.null(refSeqID)){
    prot = prot[refseq.ID == refSeqID]
    if(nrow(prot) == 0){
      stop(paste0(refSeqID, " not found!"))
    }
  } else if(!is.null(proteinID)){
    prot = prot[protein.ID == proteinID]
    if(nrow(prot) == 0){
      stop(paste0(refSeqID, " not found!"))
    }
  } else{
    txs = unique(prot$refseq.ID)
    if(length(txs) > 1){
      message(paste(length(txs), ' transcripts available. Use arguments refSeqID or proteinID to manually specify tx name.', sep = ''))
      print(prot[!duplicated(protein.ID),.(HGNC, refseq.ID, protein.ID, aa.length)])
      prot = prot[which(prot$aa.length == max(prot$aa.length)),]
      if(length(unique(prot$refseq.ID)) > 1){
        prot = prot[which(prot$refseq.ID == unique(prot[,refseq.ID])[1]),]
        message(paste('Using longer transcript', unique(prot[,refseq.ID])[1], 'for now.', sep=' '))
      } else{
        message(paste('Using longer transcript', unique(prot[,refseq.ID])[1], 'for now.', sep=' '))
      }
    }
  }
  
  #Legth of protein
  len = as.numeric(max(prot$aa.length, na.rm = TRUE))
  #Remove NA's
  #prot = prot[!is.na(Label)]
  prot = prot[,domain_lenght := End - Start][order(domain_lenght, decreasing = TRUE)][,domain_lenght := NULL]
  
  #hard coded colors for variant classification if user doesnt provide any
  
  sampleSize = as.numeric(maf@summary[ID %in% 'Samples', summary])
  mutRate = round(getGeneSummary(x = maf)[Hugo_Symbol %in% geneID, MutatedSamples]/sampleSize*100, digits = 2)
  cbioSubTitle = geneID
  if(showMutationRate){
    cbioSubTitle = substitute(paste(italic(cbioSubTitle), " : [Somatic Mutation Rate: ", mutRate, "%]"))
  }
  
  if(cBioPortal){
    vc = c("Nonstop_Mutation", "Frame_Shift_Del", "Deletion","Missense_Mutation","Frame_Shift_Del_Ins",
           "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Deletion")
    vc.cbio = c("Truncating", "Truncating", "Missense", "Truncating", "Truncating", "Truncating",
                "In-frame", "In-frame")
    names(vc.cbio) = vc
    col = grDevices::adjustcolor(col = c("black", "#33A02C", "brown"), alpha.f = 0.7)
    col = c('Truncating' = col[1], 'Missense' = col[2], 'In-frame' = col[3])
  }else{
    if(is.null(colors)){
      col = get_vcColors(alpha = 0.7, named = TRUE)
    }else{
      col = colors
    }
  }
  
  #prot.dat = prot.dat[Variant_Classification != 'Splice_Site']
  #Remove 'p.'
  prot.spl = strsplit(x = as.character(prot.dat$AAChange_), split = '.', fixed = TRUE)
  prot.conv = sapply(sapply(prot.spl, function(x) x[length(x)]), '[', 1)
  
  prot.dat[,conv := prot.conv]
  #If conversions are in HGVSp_long (default HGVSp) format, we will remove strings Ter followed by anything (e.g; p.Asn1986GlnfsTer13)
  pos = gsub(pattern = 'Ter.*', replacement = '',x = prot.dat$conv)
  
  #Following parsing takes care of most of HGVSp_short and HGVSp_long format
  pos = gsub(pattern = '[[:alpha:]]', replacement = '', x = pos)
  pos = gsub(pattern = '\\*$', replacement = '', x = pos) #Remove * if nonsense mutation ends with *
  pos = gsub(pattern = '^\\*', replacement = '', x = pos) #Remove * if nonsense mutation starts with *
  pos = gsub(pattern = '\\*.*', replacement = '', x = pos) #Remove * followed by position e.g, p.C229Lfs*18
  
  
  #pos = as.numeric(sapply(strsplit(x = pos, split = '_', fixed = TRUE), '[[', 1))
  pos = as.numeric(sapply(X = strsplit(x = pos, split = '_', fixed = TRUE), FUN = function(x) x[1]))
  prot.dat[,pos := abs(pos)]
  
  if(nrow( prot.dat[is.na(pos)]) > 0){
    message(paste('Removed', nrow( prot.dat[is.na(prot.dat$pos),]), 'mutations for which AA position was not available', sep = ' '))
    #print(prot.dat[is.na(pos)])
    prot.dat = prot.dat[!is.na(pos)]
  }
  
  prot.snp.sumamry = prot.dat[,.N, .(Variant_Classification, conv, pos)]
  colnames(prot.snp.sumamry)[ncol(prot.snp.sumamry)] = 'count'
  maxCount = max(prot.snp.sumamry$count, na.rm = TRUE)
  
  prot.snp.sumamry = prot.snp.sumamry[order(pos),]
  #prot.snp.sumamry$distance = c(0,diff(prot.snp.sumamry$pos))
  
  if(cBioPortal){
    prot.snp.sumamry$Variant_Classification = vc.cbio[as.character(prot.snp.sumamry$Variant_Classification)]
  }
  
  if(maxCount <= 5){
    prot.snp.sumamry$count2 = 1+prot.snp.sumamry$count
    lim.pos = 2:6
    lim.lab = 1:5
  }else{
    prot.snp.sumamry$count2 = 1+(prot.snp.sumamry$count * (5/max(prot.snp.sumamry$count)))
    lim.pos = prot.snp.sumamry[!duplicated(count2), count2]
    lim.lab = prot.snp.sumamry[!duplicated(count2), count]
  }
  
  if(length(lim.pos) > 6){
    lim.dat = data.table::data.table(pos = lim.pos, lab = lim.lab)
    lim.dat[,posRounded := round(pos)]
    lim.dat = lim.dat[!duplicated(posRounded)]
    lim.pos = lim.dat[,pos]
    lim.lab = lim.dat[,lab]
  }
  
  if(!defaultYaxis){
    lim.pos = c(min(lim.pos), max(lim.pos))
    lim.lab = c(min(lim.lab), max(lim.lab))
  }
  
  clusterSize = 10 #Change this later as an argument to user.
  if(repel){
    prot.snp.sumamry = repelPoints(dat = prot.snp.sumamry, protLen = len, clustSize = clusterSize)
  }else{
    prot.snp.sumamry$pos2 = prot.snp.sumamry$pos
  }
  
  xlimPos = pretty(0:max(prot$aa.length, na.rm = TRUE))
  xlimPos[length(xlimPos)] = max(prot$aa.length)
  
  # if(xlimPos[length(xlimPos)] - xlimPos[length(xlimPos)-1] <= 10){
  #   xlimPos = xlimPos[-(length(xlimPos)-1)]
  # }
  #xlimPos[length(xlimPos)] = max(as.numeric(prot$aa.length), na.rm = TRUE)
  
  #-----------------------------------
  #If user asks to label points, use ggrepel to label.
  if(!is.null(labelPos)){
    prot.snp.sumamry = data.table::data.table(prot.snp.sumamry)
    
    if(length(labelPos) == 1){
      if(labelPos != 'all'){
        prot.snp.sumamry$labThis = ifelse(test = prot.snp.sumamry$pos %in% labelPos, yes = 'yes', no = 'no')
        labDat = prot.snp.sumamry[labThis %in% 'yes']
      }else{
        labDat = prot.snp.sumamry
      }
    }else{
      prot.snp.sumamry$labThis = ifelse(test = prot.snp.sumamry$pos %in% labelPos, yes = 'yes', no = 'no')
      labDat = prot.snp.sumamry[labThis %in% 'yes']
    }
    
    if(nrow(labDat) == 0){
      message(paste0("Position ",labelPos, " doesn't seem to be mutated. Here are the mutated foci."))
      print(prot.snp.sumamry[,.(pos, conv, count, Variant_Classification)][order(pos)])
      stop()
      #return(prot.snp.sumamry[,.(mutations = sum(count)), pos][order(mutations, decreasing = TRUE)])
    }
    
    
    if(collapsePosLabel){
      uniquePos = unique(labDat[,pos2])
      labDatCollapsed = data.table::data.table()
      for(i in 1:length(uniquePos)){
        uniqueDat = labDat[pos2 %in% uniquePos[i]]
        if(nrow(uniqueDat) > 1){
          maxDat = max(uniqueDat[,count2])
          maxPos = unique(uniqueDat[,pos2])
          toLabel = uniqueDat[,conv]
          toLabel = paste(toLabel[1],paste(gsub(pattern = '^[A-z]*[[:digit:]]*', replacement = '', x = toLabel[2:length(toLabel)]), collapse = '/'), sep = '/')
          labDatCollapsed = rbind(labDatCollapsed, data.table::data.table(pos2 = maxPos, count2 = maxDat, conv = toLabel))
        }else{
          labDatCollapsed = rbind(labDatCollapsed, data.table::data.table(pos2 = uniqueDat[,pos2], count2 = uniqueDat[,count2], conv = uniqueDat[,conv]))
        }
      }
      labDat = labDatCollapsed
    }
  }
  
  
  #-----------------------------------
  #Base
  domains = unique(prot[,Label])
  domain_cols = get_domain_cols()
  
  if(length(domains) > length(domain_cols)){
    domain_cols = sample(colours(), size = length(domains), replace = FALSE)
  }
  
  domain_cols = domain_cols[1:length(domains)]
  domain_cols = grDevices::adjustcolor(col = domain_cols, alpha.f = domainAlpha)
  names(domain_cols) = domains
  
  col = col[unique(as.character(prot.snp.sumamry[,Variant_Classification]))]
  
  if(showLegend){
    lo = matrix(data = c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
    graphics::layout(mat = lo, heights = c(4, 1.25))
    par(mar = c(1, 2.5, 2, 1))
  }else{
    par(mar = c(2.5, 2.5, 2, 1))
  }
  
  
  plot(0, 0, pch = NA, ylim = c(0, 6.5), xlim = c(0, len), axes = FALSE, xlab = NA, ylab = NA)
  rect(xleft = 0, ybottom = 0.2, xright = len, ytop = 0.8, col = "#95a5a6", border = bgBorderCol)
  axis(side = 1, at = xlimPos, labels = xlimPos, lwd = 1.2, font = 1,
       cex.axis = axisTextSize[1], line = -0.4)
  axis(side = 2, at = lim.pos, labels = lim.lab, lwd = 1.2, font = 1, las = 2,
       cex.axis = axisTextSize[2])
  #mtext(text = "# Mutations", side = 2, line = 1.5, font = 1)
  segments(x0 = prot.snp.sumamry[,pos2], y0 = 0.8, x1 = prot.snp.sumamry[,pos2], y1 = prot.snp.sumamry[,count2-0.03], lwd = 1.2, col = "gray70")
  point_cols = col[as.character(prot.snp.sumamry$Variant_Classification)]
  points(x = prot.snp.sumamry[,pos2], y = prot.snp.sumamry[,count2], col = point_cols, pch = 16, cex = pointSize)
  
  prot[, domainCol := domain_cols[prot[, Label]]]
  if(roundedRect){
    if(requireNamespace("berryFunctions", quietly = TRUE)){
      for(i in 1:nrow(prot)){
        berryFunctions::roundedRect(xleft = prot[i,Start], ybottom = 0.1, xright = prot[i,End], ytop = 0.9, col = prot[i, domainCol], border = domainBorderCol, rounding = 0.08)
      }
    }else{
      #warning("Package berryFunctions needed for roundedRect to work. Please install it and try again.")
      rect(xleft = prot[,Start], ybottom = 0.1, xright = prot[,End], ytop = 0.9, col = prot[,domainCol], border = domainBorderCol)
    }
  }else{
    rect(xleft = prot[,Start], ybottom = 0.1, xright = prot[,End], ytop = 0.9, col = prot[,domainCol], border = domainBorderCol)
  }
  
  
  title(main = cbioSubTitle, adj = 0, font.main = 2, cex.main = titleSize[1], line = 0.8)
  title(main = unique(prot[,refseq.ID]), adj = 0, font.main = 1, line = -0.5, cex.main = titleSize[2])
  
  if(showDomainLabel){
    if(labelOnlyUniqueDoamins){
      prot = prot[!duplicated(Label)]
    }
    prot$pos = rowMeans(x = prot[,.(Start, End)], na.rm = FALSE)
    text(y = 0.5, x = prot$pos, labels = prot$Label, font = 3, cex = domainLabelSize)
  }
  
  if(!is.null(labelPos)){
    #prot.snp.sumamry = repelPoints(dat = prot.snp.sumamry, protLen = len, clustSize = 5)
    text(x = labDat[,pos2], y = labDat[,count2+0.45], labels = labDat[,conv],
         font = 1, srt = labPosAngle, cex = labPosSize, adj = 0.1)
  }
  
  if(showLegend){
    par(mar = c(0, 0.5, 1, 0), xpd = TRUE)
    
    plot(NULL,ylab='',xlab='', xlim=0:1, ylim=0:1, axes = FALSE)
    lep = legend("topleft", legend = names(col), col = col,  bty = "n", border=NA,
                 xpd = TRUE, text.font = 1, pch = 16, xjust = 0, yjust = 0,
                 cex = legendTxtSize, y.intersp = 1.5, x.intersp = 1,
                 pt.cex = 1.2 * legendTxtSize, ncol = ceiling(length(col)/4))
    
    x_axp = 0+lep$rect$w
    
    if(!showDomainLabel){
      if(length(domain_cols) <= 4){
        n_col = 1
      }else{
        n_col = (length(domain_cols) %/% 4)+1
      }
      
      lep = legend(x = x_axp, y = 1, legend = names(domain_cols),
                   col = domain_cols, border = NA,
                   ncol= n_col, pch = 15, xpd = TRUE, xjust = 0, bty = "n",
                   cex = legendTxtSize, title = "Domains",
                   title.adj = 0, pt.cex = 1.2 * legendTxtSize)
    }
  }
  
  if(printCount){
    print(prot.snp.sumamry[,.(pos, conv, count, Variant_Classification)][order(pos)])
    #print(prot.snp.sumamry[,.(mutations = sum(count)), pos][order(mutations, decreasing = TRUE)])
  }
}

get_domain_cols = function(){
  c("#f3a683", "#f7d794", "#778beb", "#e77f67", "#cf6a87", "#f19066",
    "#f5cd79", "#546de5", "#e15f41", "#c44569", "#786fa6", "#f8a5c2",
    "#63cdda", "#ea8685", "#596275", "#574b90", "#f78fb3", "#3dc1d3",
    "#e66767", "#303952")
}

plotmafSummary_custom = function(maf, rmOutlier = TRUE, dashboard = TRUE, titvRaw = TRUE, log_scale = FALSE,
                                 addStat = NULL, showBarcodes = FALSE, fs = 1,
                                 textSize = 0.8, color = NULL, titleSize = c(1, 0.8), titvColor = NULL, top = 10){
  
  
  addStat.opts = c('mean', 'median')
  if(!is.null(addStat)){
    if(length(addStat) > 1){
      stop('addStat can only be either mean or median.')
    }
    
    if(! addStat %in% addStat.opts){
      stop('addStat can only be either mean or median.')
    }
  }
  
  
  if(dashboard){
    #Plot in dashboard style
    pie = FALSE
    dashboard(maf = maf, color = color, rmOutlier = TRUE, log_conv = log_scale,
              titv.color = titvColor, fontSize = fs, titleSize = titleSize, sfs = statFontSize,
              n = top, donut = pie, rawcount = titvRaw, stat = addStat, barcodes = showBarcodes, barcodeSize = textSize)
    
  }else{
    
    if(is.null(color)){
      #hard coded color scheme if user doesnt provide any
      col = get_vcColors()
    }else{
      col = color
    }
    
    vcs = getSampleSummary(maf)
    vcs = vcs[,colnames(vcs)[!colnames(x = vcs) %in% c('total', 'Amp', 'Del', 'CNV_total')], with = FALSE]
    
    vcs = vcs[,c(1,order(colSums(x = vcs[,2:(ncol(vcs)), with =FALSE]), decreasing = TRUE)+1), with =FALSE] #order based on most event
    vcs.m = data.table::melt(data = vcs, id = 'Tumor_Sample_Barcode')
    colnames(vcs.m) = c('Tumor_Sample_Barcode', 'Variant_Classification', 'N')
    
    data.table::setDF(vcs)
    rownames(x = vcs) = vcs$Tumor_Sample_Barcode
    vcs = vcs[,-1]
    vcs = t(vcs)
    
    #--------------------------- variant per sample plot -----------------
    
    graphics::layout(mat = matrix(c(1, 1, 2, 2, 3, 3), nrow = 3, byrow = TRUE), heights = c(4, 4, 3))
    if(showBarcodes){
      par(mar = c(7, 2, 3, 1))
      b = barplot(vcs, col = col[rownames(vcs)], border = NA, axes = FALSE, names.arg =  rep("", ncol(vcs)))
      axis(side = 1, at = b, labels = colnames(vcs), font = 2, cex.axis = textSize, las = 2, lwd = 1)
      
    }else{
      par(mar = c(2, 2, 3, 1))
      b = barplot(vcs, col = col[rownames(vcs)], border = NA, axes = FALSE, names.arg =  rep("", ncol(vcs)))
    }
    
    axis(side = 2, at = as.integer(seq(0, max(colSums(vcs)), length.out = 4)), lwd = 2, font = 2, las = 2,
         line = -0.3, hadj = 0.6, cex.axis = fs)
    title(main = "Variants per sample", adj = 0, cex.main = titleSize[1], font = 2, line = 2)
    
    if(!is.null(addStat)){
      if(addStat == 'mean'){
        med.line = round(maf@summary[nrow(maf@summary),Mean], 2)
        df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Mean: ', med.line, sep='')))
      }else if(addStat == 'median'){
        med.line = round(maf@summary[nrow(maf@summary),Median], 2)
        df = data.frame(y = c(med.line), x = as.integer(0.8*nrow(getSampleSummary(maf))), label = c(paste('Median: ', med.line, sep='')))
      }
    }else{
      med.line = round(max(maf@summary[,Median], na.rm = TRUE), 2)
    }
    
    title(main = paste0("Median: ", med.line), adj = 0, cex.main = titleSize[2], font = 2, line = 1)
    lines(x = c(1, b[length(b)]), y = c(med.line, med.line), col = "maroon", lwd = 2, lty = 2)
    
    #--------------------------- vc summary plot -----------------
    par(mar = c(2, 2, 2, 1))
    boxH = vcs.m[,boxplot.stats(N)$stat[5], by = .(Variant_Classification)]
    colnames(boxH)[ncol(boxH)] = 'boxStat'
    b = boxplot(N ~ Variant_Classification, data = vcs.m, col = col[levels(vcs.m$Variant_Classification)],
                axes = FALSE, outline = FALSE, lwd = 1, border = grDevices::adjustcolor(col = "black", alpha.f = 0.6))
    axis(side = 2, at = as.integer(seq(0, max(boxH[,boxStat], na.rm = TRUE), length.out = 4)),
         lwd = 2, font = 2, cex.axis = fs, las = 2)
    title(main = "Variant Classification summary", adj = 0, cex.main = fs, font = 2, line = 1)
    
    plot.new()
    par(mar = c(4, 2.5, 0.5, 2))
    legend(x = "top", legend = names(col[levels(vcs.m$Variant_Classification)]),
           fill = col[levels(vcs.m$Variant_Classification)],
           bty = "n", ncol = 3, cex = fs)
  }
}

