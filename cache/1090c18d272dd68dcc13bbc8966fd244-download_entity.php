<!DOCTYPE html>
<html>
<head>
<link rel="stylesheet" type="text/css" href="CSS/signor_style.css">
<link rel="shortcut icon" href="favicon.ico" />
<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.0/jquery.min.js"></script>

<title>Downloads</title>

    
    
<meta charset="UTF-8">
</head>

<body>
<script>
  function causalTabx() {
    var element=document.getElementById("human");
    var causal=document.getElementById("causalTab");
    if (element && element.checked === true) {
      if (causal){
        causal.disabled=false;
      }
      
    }
    else{
      if (causal){
        causal.disabled=true;
      }
    }
  }
  $(document).ready(function(){
      causalTabx();
    });
  
</script>
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-85188337-1', 'auto');
  ga('send', 'pageview');
  

</script>
 
  <div class=whole_wrapper style='background: rgb(241, 241, 240)'>
<header>
   <div style="text-align: center;" class="logo_area">
      <img src="img/signor_logo.png" style="display: inline-block;" class="head_logo" alt="Signor logo">
      <img src="img/signor_version.png" class="head_logo" style="display: inline-block;">
      <h1>The
         <span class="header_letters">SIG</span>naling
         <span class="header_letters">N</span>etwork
         <span class="header_letters">O</span>pen
         <span class="header_letters">R</span>esource
      </h1>
   </div>
   <ul class="menu">
      <li class="menu_list"><a class="active" href="http://signor.uniroma2.it"><b>HOME</b></a></li>
      <li class="menu_list"><a href="user_guide.php"><b>USER GUIDE</b></a></li>
      <li class="menu_list"><a href="statistics.php"><b>STATISTICS</b></a></li>
      <li class="menu_list"><a href="curation.php"><b>CURATION</b></a></li>
      <li class="menu_list"><a href="contact.php"><b>CONTACT</b></a></li>
      <li class="menu_list"><a href="downloads.php"><b>DOWNLOADS</b></a></li>
      <li class="menu_list"><a href="APIs.php"><b>APIs</b></a></li>
      
      <li class="menu_list">
         <form class=search_entity style="margin-bottom:0px;" action="search_result.php" method="post" enctype="multipart/form-data">
            <span style="color:white;">
            <input type="text" style="width:200px;font-size:14px;" name="entity" placeholder="Type Entity Here..." required>
            <input class=submit_button_top type="submit" value="SEARCH">
                        </span>
         </form>
          
      </li>
   </ul>
</header>
<div class=content style="margin: auto; max-width: 1300px">
<div class='menu_down'>  
  <div class=left_downloads style='background: white;width:90%; min-width: 305px;'>
    <h4 class=div_title_home>Available Downloads</h4>
    <div class=downloads_list>
      
     <b>
    <table style='margin:auto; width:90%; border-collapse: collapse;'>
      
      <tr style='line-height: 2;border-bottom: green 2px solid;'>
        <td style='text-align: left'>
          Whole database data
        </td>
        <td>
          <a class='download_jump' href='#all_download'>view</a>
        </td>
      </tr>
      <tr style='line-height: 2;border-bottom: green 2px solid;'>
        <td style='text-align: left'>
          Programmatic Access
        </td>
        <td>
          <a class='download_jump' href='/APIs.php' target="_blank">view</a>
        </td>
      </tr>
       <tr style='line-height: 3;border-bottom: green 2px solid;'>
        <td style='text-align: left'>
          SIGNOR entity data
        </td>
        <td>
          <a class='download_jump' href='#signor_entity'>view</a>
        </td>
      </tr>
      <tr style='line-height: 3;border-bottom: green 2px solid;'>
        <td style='text-align: left'>
          Pathway-related data
        </td>
        <td>
          <a class='download_jump' href='#path_download'>view</a>
        </td>
      </tr>
      <tr style='line-height: 3;border-bottom: green 2px solid;'>
        <td style='text-align: left'>
          Phosphorylation data
        </td>
        <td>
          <a class='download_jump' href='#phospho_data'>view</a>
        </td>
      </tr>
      <tr style='line-height: 3;border-bottom: green 2px solid;'>
        <td style='text-align: left'>
          Download by relation id
        </td>
        <td>
          <a class='download_jump' href='#relation'>view</a>
        </td>
      </tr>
      <tr style='line-height: 3;'>
        <td style='text-align: left'>
          Download by entity or list of entities
        </td>
        <td>
          <a class='download_jump' href='#entities'>view</a>
        </td>
      </tr>
      
      
    </table>
   </b> 
    </div>
    
  </div>
  <div style='display: inline-block'>
      <span class=indicator style='display: none; color:red;'>&#9658;</span>
  </div>
</div>
  <div class=right_downloads>
    <div class=downloads>
    
    <h4 class=div_title_home id='all_download'>Download All Data</h4>
    
    <form action="download_entity.php" method="post" enctype="multipart/form-data">

    
    <b>Select Organism:</b><br><br>
     <table class=organism_table>
      <tr>
        <td>
          <input id="human" onchange="causalTabx()" type="radio" name="organism" value="human" checked>
        </td>
        <td style='text-align: center'>
          <img width='15px'src='img/man.png'>
        </td>
        <td>
          Homo sapiens
        </td>
      </tr>
      <tr>
        <td>
          <input type="radio" onchange="causalTabx()" name="organism" value="mouse">
        </td>
        <td style='text-align: center'>
          <img width='40px'src='img/mouse.png'>
        </td>
        <td>
          Mus musculus 
        </td>
      </tr>
      <tr>
        <td>
          <input type='radio' onchange="causalTabx();" name='organism' value='rat'>
        </td>
        <td style='text-align: center'>
          <img width='30px'src='img/rat.png'>
        </td>
        <td>
          Rattus norvegicus
        </td>
      </tr>
    </table>
   
     <h4>Select a format below to download all data from the database:</h4>
    <input type="radio" name="format" value="csv">file tsv <br>
    <input type="radio" name="format" value="Excel5">file xls (slower)<br>
    <input id="causalTab" disabled=true type="radio" name="format" value="causalTab">causalTab (BETA) - <a onClick="window.open('scripts/causalTabInfo.php','causalTab Info','resizable,height=450px,width=450px'); return false;" style="font-size:small;color: #52a677;text-decoration: underline;cursor:pointer">more info</a><br>
    <br>
    <input type="submit" onclick="ga('send', 'event', 'Downloads', 'Click', 'Data downloaded');" value="Download" name="submit">
    
    
  </form>
    </div>
  <div class=downloads id='signor_entity'>
    <h4 class=div_title_home>Download SIGNOR entity data</h4>
    <form action="download_complexes.php" method="post" enctype="multipart/form-data">
    <b>SIGNOR Complexes:</b><br><br><input type="submit" value="Download complex data" name="submit">
    </form>
    <form action="download_complexes.php" method="post" enctype="multipart/form-data">
    <b>SIGNOR Protein Families:</b><br><br><input type="submit" value="Download protein family data" name="submit">
    </form>
    <form action="download_signor_def.php" method="post" enctype="multipart/form-data">
      <b>SIGNOR Controlled Vocabulary</b><br><br>
      <table>
        <tr>
          <td>Entity Types:</td>
          <td><input type="submit" value="Download Entity-Type CV" name="submit"></td>
        </tr>
        <tr>
          <td>
            Modification Types:
          </td>
          <td>
            <input type="submit" value="Download Modifications CV" name="submit">
          </td>
        </tr>
        <tr>
          <td>Mechanisms Types:</td>
          <td><input type="submit" value="Download Mechansims CV" name="submit"></td>
        </tr>
      </table>

    </form>
  </div>
    <div class=downloads_pathway id='path_download'>
    <h4 class=div_title_home>Download Pathway Data</h4>
 
    <form action="download_entity.php" method="post" enctype="multipart/form-data">
          Select a pathway from the list below to download its data:
          <br><br><br>
         <select required name=pathway_list><option value=''>Select option below</option><option value='SIGNOR-AML'>Acute Myeloid Leukemia</option><option value='SIGNOR-AC'>Adipogenesis</option><option value='SIGNOR-Adipogenesis: BMP'>Adipogenesis: BMP</option><option value='SIGNOR-AD'>Alzheimer Disease</option><option value='SIGNOR-AMPK'>AMPK Signaling</option><option value='SIGNOR-Autophagy'>Autophagy</option><option value='SIGNOR-G1-S_trans'>Cell cycle: G1/S phase transition</option><option value='SIGNOR-G2-M_trans'>Cell cycle: G2/M phase transition</option><option value='SIGNOR-CRC'>Colorectal Carcinoma</option><option value='SIGNOR-DR'>Death Receptor Signaling</option><option value='SIGNOR-EcmSynthesis'>ECM: synthesis</option><option value='SIGNOR-EGF'>EGFR</option><option value='SIGNOR-EosCCL11'>Eosinophil: CCL11</option><option value='SIGNOR-EosCCR1'>Eosinophil: CCR1</option><option value='SIGNOR-EosIL5'>Eosinophil: IL5</option><option value='SIGNOR-FapATC'>FAP: adipogenic transcriptional cascade</option><option value='SIGNOR-FapBMP'>FAP: BMP</option><option value='SIGNOR-FapGCC'>FAP: glucocorticoids</option><option value='SIGNOR-FapIBMXR'>FAP: IBMX/rosiglitazone</option><option value='SIGNOR-FapIL4'>FAP: IL4</option><option value='SIGNOR-FapINS'>FAP: insulin-mediated adipogenesis</option><option value='SIGNOR-FapNOTCH'>FAP: Notch</option><option value='SIGNOR-FapPDGFRA'>FAP: PDGFR alpha</option><option value='SIGNOR-FapSHH'>FAP: SHH</option><option value='SIGNOR-FapSIR'>FAP: SIRT1</option><option value='SIGNOR-FapTGFB'>FAP: TGF beta</option><option value='SIGNOR-FapWNT'>FAP: Wnt</option><option value='SIGNOR-FibroFN1TR'>Fibroblast: FN1 transcriptional regulation</option><option value='SIGNOR-FSGS'>FSGS</option><option value='SIGNOR-GBM'>Glioblastoma Multiforme</option><option value='SIGNOR-GCR'>Glucocorticoid Receptor</option><option value='SIGNOR-Hedgehog'>Hedgehog</option><option value='SIGNOR-HPP'>Hippo Signaling</option><option value='SIGNOR-Myogenesis'>IGF and Myogenesis</option><option value='SIGNOR-IL1R'>IL1 Receptor</option><option value='SIGNOR-IL6'>IL6 </option><option value='SIGNOR-IOA'>Inhibition of Apoptosis</option><option value='SIGNOR-INSR'>Insulin Receptor</option><option value='SIGNOR-LBC'>Luminal Breast Cancer</option><option value='SIGNOR-MacroCCL2TR'>Macrophage: CCL2 transcriptional regulation</option><option value='SIGNOR-M1M2'>Macrophage differentiation</option><option value='SIGNOR-MacroIFNG'>Macrophage: IFN gamma</option><option value='SIGNOR-MacroIL10'>Macrophage: IL10</option><option value='SIGNOR-MacroIL10TR'>Macrophage: IL10 transcriptional regulation</option><option value='SIGNOR-MacroIL1'>Macrophage: IL1 beta</option><option value='SIGNOR-MacroIL1BTR'>Macrophage: IL1 beta transcriptional regulation</option><option value='SIGNOR-MacroIL4'>Macrophage: IL4</option><option value='SIGNOR-MacroIL6'>Macrophage: IL6</option><option value='SIGNOR-MacroIL6TR'>Macrophage: IL6 transcriptional regulation</option><option value='SIGNOR-MT'>Macrophage: TGF beta</option><option value='SIGNOR-MacroTGFBTR'>Macrophage: TGF beta transcriptional regulation</option><option value='SIGNOR-MacroTNFA'>Macrophage: TNF alpha</option><option value='SIGNOR-MacroTNFATR'>Macrophage: TNF alpha transcriptional regulation</option><option value='SIGNOR-MM'>Malignant Melanoma</option><option value='SIGNOR-MastFCERI'>Mast cell: FCERI</option><option value='SIGNOR-MastSCF'>Mast cell: SCF</option><option value='SIGNOR-MCAPO'>Mitochondrial Control of Apoptosis</option><option value='SIGNOR-MonoIL10TR'>Monocyte: IL10 transcriptional regulation</option><option value='SIGNOR-MS'>MTOR Signaling</option><option value='SIGNOR-IL6HYPER'>Myoblast: IL6 and hypertrophy</option><option value='SIGNOR-MMDD'>Myogenesis modulation: DNA damage</option><option value='SIGNOR-Myostatin'>Myostatin</option><option value='SIGNOR-MyoATR'>Myotubes: atrophy transcriptional regulation</option><option value='SIGNOR-MyoIL6ATR'>Myotubes: IL6 and atrophy</option><option value='SIGNOR-NeutroAM'>Neutrophils: adhesion and migration</option><option value='SIGNOR-NeutroGPCR'>Neutrophils: G-protein-coupled receptor</option><option value='SIGNOR-NFKBC'>NF-KB Canonical</option><option value='SIGNOR-NFKBNC'>NF-KB Non Canonical</option><option value='SIGNOR-NS'>Noonan syndrome</option><option value='SIGNOR-NOTCH'>NOTCH Signaling</option><option value='SIGNOR-NOTCH_Myogenesis'>NOTCH Signaling and Myogenesis</option><option value='SIGNOR-OstCAL'>Osteoblast: osteocalcin transcriptional regulation</option><option value='SIGNOR-P38'>P38 Signaling</option><option value='SIGNOR-P38_Myogenesis'>P38 Signaling and Myogenesis</option><option value='SIGNOR-PD'>Parkinson disease</option><option value='SIGNOR-PI3K/AKT'>PI3K/AKT</option><option value='SIGNOR-PC'>Prostate Cancer</option><option value='SIGNOR-RMS'>Rhabdomyosarcoma</option><option value='SIGNOR-SAPK-JNK'>SAPK/JNK Signaling</option><option value='SIGNOR-ScATR'>Satellite: atrophy</option><option value='SIGNOR-ScCA'>Satellite: calcium-mediated differentiation</option><option value='SIGNOR-ScCCR2'>Satellite: CCR2</option><option value='SIGNOR-ScFSTTR'>Satellite: FST transcriptional regulation</option><option value='SIGNOR-ScIFNG'>Satellite: IFN gamma</option><option value='SIGNOR-ScIGF1'>Satellite: IGF1</option><option value='SIGNOR-ScIL1SEC'>Satellite: IL1 and secretome</option><option value='SIGNOR-ScIL6'>Satellite: IL6</option><option value='SIGNOR-ScMTC'>Satellite: myogenic transcriptional cascade</option><option value='SIGNOR-ScNOTCH'>Satellite: NOTCH</option><option value='SIGNOR-SPAM'>Satellite: prostaglandins and myogenesis</option><option value='SIGNOR-SateTNFA'>Satellite: TNFA</option><option value='SIGNOR-ScWNT'>Satellite: WNT</option><option value='SIGNOR-TCA'>T cell activation</option><option value='SIGNOR-TGFb'>TGFbeta Signaling</option><option value='SIGNOR-TC'>Thyroid cancer</option><option value='SIGNOR-TLymphATP'>T-Lymphocyte: ATP</option><option value='SIGNOR-TLymphIFNG'>T-Lymphocyte: IFN gamma</option><option value='SIGNOR-TLymphIFNGTR'>T-Lymphocyte: IFN gamma transcriptional regulation</option><option value='SIGNOR-TLymphIL10'>T-Lymphocyte: IL10</option><option value='SIGNOR-TLymphIL10TR'>T-Lymphocyte:IL10 transcriptional regulation</option><option value='SIGNOR-TLymphIL4'>T-Lymphocyte: IL4</option><option value='SIGNOR-TLymphIL4TC'>T-Lymphocyte: IL4 transcriptional regulation</option><option value='SIGNOR-TLymphIL6'>T-Lymphocyte: IL6</option><option value='SIGNOR-TLITR'>T-lymphocytes: IL10 transcriptional regulation</option><option value='SIGNOR-TLymphTGFB'>T-Lymphocyte: TGF beta</option><option value='SIGNOR-TLymphTGFBTR'>T-lymphocyte: TGF beta transcriptional regulation</option><option value='SIGNOR-TA'>TNF alpha</option><option value='SIGNOR-Tnfa_Notch_crosstalk'>Tnfa_Notch_crosstalk</option><option value='SIGNOR-TLR'>Tol like receptors</option><option value='SIGNOR-WNT'>WNT Signaling</option><option value='SIGNOR-WNT_Myogenesis'>WNT Signaling and Myogenesis</option></select>    <input type="submit" value="Download" name="submit">
    <input type="submit"  value="Download SBML (BETA)" name="sbml">
    
  </form>
    </div>
    <div class=downloads id=phospho_data>
    <h4 class=div_title_home>Phosphorylation Data</h4>
    <form action="download_entity.php" method="post" enctype="multipart/form-data">

    <b>Select Organism:</b><br><br>
    <table class=organism_table>
      <tr>
        <td>
          <input type="radio" name="organism" value="human" checked>
        </td>
        <td style='text-align: center'>
          <img width='15px'src='img/man.png'>
        </td>
        <td>
          Homo sapiens
        </td>
      </tr>
      <tr>
        <td>
          <input type="radio" name="organism" value="mouse">
        </td>
        <td style='text-align: center'>
          <img width='40px'src='img/mouse.png'> 
        </td>
        <td>
          Mus musculus
        </td>
      </tr>
      <tr>
        <td>
          <input type='radio' name='organism' value='rat'>
        </td>
        <td style='text-align: center'>
          <img width='30px'src='img/rat.png'>
        </td>
        <td>
          Rattus norvegicus
        </td>
      </tr>
    </table>
     
   
     <h4>Select a format below to download phosphorylation data from the database:</h4>
     
    <input type="radio" name="format" value="csv">file tsv <br>
    <input type="radio" name="format" value="Excel5">file xls (slower)<br>
    
    <br>
     <input type="hidden" name="phosphorylation" value="phosphorylation">
    <input type="submit" value="Download" name="submit">
    </form>
  </div>
    <div class=downloads_relation id='relation'>
    <h4 class=div_title_home>Download Relations</h4>
 
    <form action="download_entity.php" method="post" enctype="multipart/form-data">
    <input style="width:200px; height:30px; font-size:medium;" type="text" name="relation_list" placeholder="Type Relation ID..." required><br><br>
    <input type="radio" name="format" value="csv" required>file tsv <br>
    <input type="radio" name="format" value="Excel5" required>file xls (slower)<br> <br>
    <input type="submit" value="Download" name="submit">
    </form>
    </div>
    

     
  
    <div class=downloads_entity id='entities'>
    <h4 class=div_title_home>Download Entity Data</h4>
    <form action="download_entity.php" method="post" enctype="multipart/form-data">
     <b>Download instructions:</b><br><br>        
            <ol>

      <li>Type entity IDs separated by a delimiter.<br><br></li>
      <li>Select 'all' to download relation data involving the list of entities provided.<br><br></li>
      <li>Select 'connect' to download data of relations <b>only</b> among the entities in the list.<br><br></li>
      <li>Choose a file format.</li>
      </ol>
    <input style="width:300px; height:30px; font-size:medium;" type="text" name="list" placeholder="Type Entity Here..." required>
    <input type="radio" name="data_type" value="all">all<input type="radio" name="data_type" value="connect">connect<br><br>
    <input type="radio" name="format" value="csv" required>file tsv <br>
    <input type="radio" name="format" value="Excel5" required>file xls (slower)<br> <br>
    <input type="submit" value="Download" name="submit">
    </form>
    </div>
    
      </div>
       

  </div>
  
<script>
  $('.download_jump').on('click', function(){
    $('.indicator').css({
      'display' : '',
      'position' : '',
      'bottom' : ''
      });
    $('.indicator').show();

  });
  $('.download_jump:last').click(function(){
    $('.indicator').hide();
    $('.indicator').css({
      'display' : 'inline-block',
      'position' : 'absolute',
      'bottom' : '15px'
      });
  });
</script>
<footer>
    © 2018 SIGNOR
</footer>
  </div>
</body><script type='text/javascript'>alert('No relations found between the selected entities!');</script>