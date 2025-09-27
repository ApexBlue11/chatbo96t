import streamlit as st
import requests
import pandas as pd
import json
import time
import re
from typing import Dict, Any, List, Optional
from dataclasses import dataclass
from io import StringIO

# Configure Streamlit page
st.set_page_config(
    page_title="Genetic Variant Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
.main-header {
    font-size: 2.5rem;
    color: #1f77b4;
    text-align: center;
    margin-bottom: 2rem;
}
.section-header {
    font-size: 1.5rem;
    color: #2c3e50;
    border-bottom: 2px solid #3498db;
    padding-bottom: 0.5rem;
    margin: 1rem 0;
}
.success-box {
    background-color: #d4edda;
    border: 1px solid #c3e6cb;
    color: #155724;
    padding: 1rem;
    border-radius: 0.5rem;
    margin: 1rem 0;
}
.error-box {
    background-color: #f8d7da;
    border: 1px solid #f5c6cb;
    color: #721c24;
    padding: 1rem;
    border-radius: 0.5rem;
    margin: 1rem 0;
}
.info-box {
    background-color: #d1ecf1;
    border: 1px solid #bee5eb;
    color: #0c5460;
    padding: 1rem;
    border-radius: 0.5rem;
    margin: 1rem 0;
}
</style>
""", unsafe_allow_html=True)

@dataclass
class QueryClassification:
    is_genomic: bool
    query_type: str
    extracted_identifier: Optional[str]

class GenomicQueryRouter:
    def __init__(self):
        self.hgvs_patterns = {
            'transcript': [
                r'\b(NM_\d+(?:\.\d+)?):c\.[A-Za-z0-9\-+*>_]+',
                r'\b(ENST\d+(?:\.\d+)?):c\.[A-Za-z0-9\-+*>_]+',
            ],
            'genomic': [
                r'\b(NC_\d+(?:\.\d+)?):g\.[A-Za-z0-9\-+*>_]+',
                r'\b(chr(?:\d+|X|Y|MT?)):g\.\d+[A-Za-z]+>[A-Za-z]+',
            ],
            'protein': [
                r'\b(NP_\d+(?:\.\d+)?):p\.[A-Za-z0-9\-+*>_()]+',
                r'\b(ENSP\d+(?:\.\d+)?):p\.[A-Za-z0-9\-+*>_()]+',
            ]
        }
        self.rsid_pattern = r'\b(rs\d+)\b'
    
    def classify_query(self, query: str) -> QueryClassification:
        query = query.strip()
        
        for variant_type, patterns in self.hgvs_patterns.items():
            for pattern in patterns:
                match = re.search(pattern, query, re.IGNORECASE)
                if match:
                    return QueryClassification(
                        is_genomic=True,
                        query_type=f'hgvs_{variant_type}',
                        extracted_identifier=match.group(0)
                    )
        
        rsid_match = re.search(self.rsid_pattern, query, re.IGNORECASE)
        if rsid_match:
            return QueryClassification(
                is_genomic=True,
                query_type='rsid',
                extracted_identifier=rsid_match.group(1)
            )
        
        return QueryClassification(
            is_genomic=False,
            query_type='general',
            extracted_identifier=None
        )

def query_clingen_allele(hgvs: str) -> Dict[str, Any]:
    """Query ClinGen Allele Registry by HGVS notation."""
    base_url = "http://reg.clinicalgenome.org/allele"
    params = {'hgvs': hgvs}
    
    with st.spinner(f"Querying ClinGen for: {hgvs}"):
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()

def parse_caid_minimal(raw_json):
    """Parse ClinGen Allele Registry JSON to extract key information."""
    result = {}

    # CAid - extract from @id URL
    result['CAid'] = raw_json.get('@id', '').split('/')[-1]

    # RSID from dbSNP external records
    dbsnp = raw_json.get('externalRecords', {}).get('dbSNP', [])
    result['rsid'] = dbsnp[0].get('rs') if dbsnp else None

    # Genomic HGVS for both GRCh38 and GRCh37
    genomic = raw_json.get('genomicAlleles', [])
    result['genomic_hgvs_grch38'] = None
    result['genomic_hgvs_grch37'] = None
    
    for g in genomic:
        hgvs_list = g.get('hgvs', [])
        ref_genome = g.get('referenceGenome', '')
        
        if 'GRCh38' in ref_genome and hgvs_list:
            result['genomic_hgvs_grch38'] = hgvs_list[0]
        elif 'GRCh37' in ref_genome and hgvs_list:
            result['genomic_hgvs_grch37'] = hgvs_list[0]

    # MyVariantInfo IDs for both builds
    myv = raw_json.get('externalRecords', {})
    result['myvariant_hg38'] = myv.get('MyVariantInfo_hg38', [{}])[0].get('id') if myv.get('MyVariantInfo_hg38') else None
    result['myvariant_hg19'] = myv.get('MyVariantInfo_hg19', [{}])[0].get('id') if myv.get('MyVariantInfo_hg19') else None

    # MANE Select transcripts (both Ensembl and RefSeq)
    result['mane_ensembl'] = None
    result['mane_refseq'] = None
    
    transcripts = raw_json.get('transcriptAlleles', [])
    for t in transcripts:
        mane = t.get('MANE', {})
        if mane and mane.get('maneStatus') == 'MANE Select':
            if 'nucleotide' in mane:
                ensembl_info = mane['nucleotide'].get('Ensembl', {})
                refseq_info = mane['nucleotide'].get('RefSeq', {})
                
                result['mane_ensembl'] = ensembl_info.get('hgvs')
                result['mane_refseq'] = refseq_info.get('hgvs')
            break

    return result

def get_variant_annotations(clingen_data):
    """Retrieve variant annotations from multiple APIs."""
    annotations = {
        'myvariant_data': {},
        'vep_data': [],
        'errors': []
    }
    
    # MyVariant.info query
    if clingen_data.get('myvariant_hg38'):
        try:
            with st.spinner("Fetching MyVariant.info data..."):
                myv_url = f"https://myvariant.info/v1/variant/{clingen_data['myvariant_hg38']}?assembly=hg38"
                myv_response = requests.get(myv_url, timeout=30)
                if myv_response.ok:
                    annotations['myvariant_data'] = myv_response.json()
                else:
                    annotations['errors'].append(f"MyVariant query failed: HTTP {myv_response.status_code}")
        except Exception as e:
            annotations['errors'].append(f"MyVariant query error: {str(e)}")
    
    # Ensembl VEP query
    if clingen_data.get('mane_ensembl'):
        try:
            with st.spinner("Fetching Ensembl VEP data..."):
                vep_input = clingen_data['mane_ensembl']
                vep_url = f"https://rest.ensembl.org/vep/human/hgvs/{vep_input}"
                vep_headers = {"Content-Type": "application/json", "Accept": "application/json"}
                vep_response = requests.get(vep_url, headers=vep_headers, timeout=30)
                if vep_response.ok:
                    annotations['vep_data'] = vep_response.json()
                else:
                    annotations['errors'].append(f"VEP query failed: HTTP {vep_response.status_code}")
        except Exception as e:
            annotations['errors'].append(f"VEP query error: {str(e)}")
    
    return annotations

def create_frequency_chart(myvariant_data):
    """Create a population frequency visualization using Streamlit's built-in charting."""
    freq_data = []
    
    # Extract frequencies from multiple possible locations in the data
    
    # 1. Check gnomAD exome data
    gnomad_exome = myvariant_data.get('gnomad_exome', {})
    if isinstance(gnomad_exome, dict) and 'af' in gnomad_exome:
        af_data = gnomad_exome['af']
        if isinstance(af_data, dict):
            for pop, freq in af_data.items():
                if isinstance(freq, (int, float)) and freq > 0:
                    clean_name = pop.replace('af_', '').upper()
                    freq_data.append({'Population': f'gnomAD Exome {clean_name}', 'Frequency': freq})
    
    # 2. Check gnomAD genome data  
    gnomad_genome = myvariant_data.get('gnomad_genome', {})
    if isinstance(gnomad_genome, dict) and 'af' in gnomad_genome:
        af_data = gnomad_genome['af']
        if isinstance(af_data, dict):
            for pop, freq in af_data.items():
                if isinstance(freq, (int, float)) and freq > 0:
                    clean_name = pop.replace('af_', '').upper()
                    freq_data.append({'Population': f'gnomAD Genome {clean_name}', 'Frequency': freq})
    
    # 3. Check dbNSFP 1000 Genomes data
    dbnsfp = myvariant_data.get('dbnsfp', {})
    if isinstance(dbnsfp, dict):
        # 1000 Genomes Project data
        kg_data = dbnsfp.get('1000gp3', {})
        if isinstance(kg_data, dict):
            # Overall frequency
            if 'af' in kg_data and isinstance(kg_data['af'], (int, float)) and kg_data['af'] > 0:
                freq_data.append({'Population': '1000G Overall', 'Frequency': kg_data['af']})
            
            # Population-specific frequencies
            for pop in ['afr', 'amr', 'eas', 'eur', 'sas']:
                if pop in kg_data and isinstance(kg_data[pop], dict):
                    pop_freq = kg_data[pop].get('af')
                    if isinstance(pop_freq, (int, float)) and pop_freq > 0:
                        freq_data.append({'Population': f'1000G {pop.upper()}', 'Frequency': pop_freq})
        
        # ExAC data
        exac_data = dbnsfp.get('exac', {})
        if isinstance(exac_data, dict):
            if 'af' in exac_data and isinstance(exac_data['af'], (int, float)) and exac_data['af'] > 0:
                freq_data.append({'Population': 'ExAC Overall', 'Frequency': exac_data['af']})
            
            for pop in ['afr', 'amr', 'eas', 'fin', 'nfe', 'sas']:
                if pop in exac_data:
                    pop_freq = exac_data[pop].get('af') if isinstance(exac_data[pop], dict) else exac_data[pop]
                    if isinstance(pop_freq, (int, float)) and pop_freq > 0:
                        freq_data.append({'Population': f'ExAC {pop.upper()}', 'Frequency': pop_freq})
    
    # 4. Check any other frequency fields in the root
    for key, value in myvariant_data.items():
        if 'af' in key.lower() and 'gnomad' in key.lower() and isinstance(value, (int, float)) and value > 0:
            clean_key = key.replace('gnomad_', '').replace('_af', '').replace('_', ' ').title()
            freq_data.append({'Population': clean_key, 'Frequency': value})
    
    if freq_data:
        # Sort by frequency for better visualization
        freq_data.sort(key=lambda x: x['Frequency'], reverse=True)
        df_freq = pd.DataFrame(freq_data)
        return df_freq
    return None

def display_clinical_significance(myvariant_data):
    """Display ClinVar clinical significance information."""
    clinvar = myvariant_data.get('clinvar', {})
    if not clinvar:
        return None
    
    clinical_info = {}
    
    # Extract key clinical information from different possible structures
    clinical_significance = (clinvar.get('clinical_significance') or 
                           clinvar.get('clnsig') or
                           clinvar.get('clinicalsignificance'))
    
    if clinical_significance:
        clinical_info['Clinical Significance'] = clinical_significance
    
    # Review status
    review_status = (clinvar.get('review_status') or 
                    clinvar.get('reviewstatus') or
                    clinvar.get('review'))
    if review_status:
        clinical_info['Review Status'] = review_status
    
    # RCV accessions - handle the array structure from your data
    rcv_data = clinvar.get('rcv', [])
    if rcv_data and isinstance(rcv_data, list):
        # Extract accession numbers and clinical significance from each RCV
        accessions = []
        significances = []
        conditions = []
        
        for rcv in rcv_data:
            if isinstance(rcv, dict):
                if rcv.get('accession'):
                    accessions.append(rcv['accession'])
                if rcv.get('clinical_significance'):
                    significances.append(rcv['clinical_significance'])
                if rcv.get('conditions', {}).get('name'):
                    conditions.append(rcv['conditions']['name'])
        
        if accessions:
            clinical_info['RCV Accessions'] = ', '.join(accessions)
        if significances:
            # Get unique significances
            unique_sigs = list(set(significances))
            clinical_info['Clinical Significance'] = ', '.join(unique_sigs)
        if conditions:
            unique_conditions = list(set(conditions))
            clinical_info['Associated Conditions'] = ', '.join(unique_conditions)
    
    # Variation ID and Allele ID
    if clinvar.get('variation_id'):
        clinical_info['Variation ID'] = clinvar['variation_id']
    if clinvar.get('allele_id'):
        clinical_info['Allele ID'] = clinvar['allele_id']
    
    # Gene information
    gene_info = clinvar.get('gene', {})
    if isinstance(gene_info, dict):
        if gene_info.get('symbol'):
            clinical_info['Gene Symbol'] = gene_info['symbol']
        if gene_info.get('id'):
            clinical_info['Gene ID'] = gene_info['id']
    
    # HGVS notations
    hgvs_info = clinvar.get('hgvs', {})
    if isinstance(hgvs_info, dict):
        if hgvs_info.get('coding'):
            clinical_info['HGVS Coding'] = hgvs_info['coding']
        if hgvs_info.get('protein'):
            clinical_info['HGVS Protein'] = hgvs_info['protein']
        if hgvs_info.get('genomic'):
            genomic = hgvs_info['genomic']
            if isinstance(genomic, list):
                clinical_info['HGVS Genomic'] = ', '.join(genomic)
            else:
                clinical_info['HGVS Genomic'] = str(genomic)
    
    return clinical_info if clinical_info else None

def display_prediction_scores(myvariant_data):
    """Display functional prediction scores."""
    dbnsfp = myvariant_data.get('dbnsfp', {})
    if not dbnsfp:
        return None
    
    predictions = {}
    
    def safe_extract_value(data, key):
        """Safely extract a value, handling lists by taking the first element."""
        if key not in data:
            return None
        value = data[key]
        if isinstance(value, list):
            return value[0] if value else None
        return value
    
    # SIFT
    if 'sift' in dbnsfp:
        sift_data = dbnsfp['sift']
        if isinstance(sift_data, dict):
            score = safe_extract_value(sift_data, 'score')
            pred = safe_extract_value(sift_data, 'pred')
            if score is not None:
                predictions['SIFT Score'] = score
            if pred is not None:
                predictions['SIFT Prediction'] = pred
    
    # Handle direct dbnsfp fields (like in your error)
    score = safe_extract_value(dbnsfp, 'sift_score')
    pred = safe_extract_value(dbnsfp, 'sift_pred')
    if score is not None:
        predictions['SIFT Score'] = score
    if pred is not None:
        predictions['SIFT Prediction'] = pred
    
    # PolyPhen
    if 'polyphen2' in dbnsfp:
        pp2_data = dbnsfp['polyphen2']
        if isinstance(pp2_data, dict):
            if 'hdiv' in pp2_data and isinstance(pp2_data['hdiv'], dict):
                score = safe_extract_value(pp2_data['hdiv'], 'score')
                pred = safe_extract_value(pp2_data['hdiv'], 'pred')
                if score is not None:
                    predictions['PolyPhen2 HDiv Score'] = score
                if pred is not None:
                    predictions['PolyPhen2 HDiv Prediction'] = pred
    
    # Handle direct polyphen fields
    pp_score = safe_extract_value(dbnsfp, 'polyphen2_hdiv_score')
    pp_pred = safe_extract_value(dbnsfp, 'polyphen2_hdiv_pred')
    if pp_score is not None:
        predictions['PolyPhen2 Score'] = pp_score
    if pp_pred is not None:
        predictions['PolyPhen2 Prediction'] = pp_pred
    
    # CADD
    if 'cadd' in dbnsfp:
        cadd_data = dbnsfp['cadd']
        if isinstance(cadd_data, dict):
            score = safe_extract_value(cadd_data, 'phred')
            if score is not None:
                predictions['CADD Phred Score'] = score
    
    # Handle direct CADD field
    cadd_score = safe_extract_value(dbnsfp, 'cadd_phred')
    if cadd_score is not None:
        predictions['CADD Phred Score'] = cadd_score
    
    return predictions if predictions else None

def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 Genetic Variant Analyzer</h1>', unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.markdown("### About")
        st.write("""
        This tool analyzes genetic variants using multiple genomic databases:
        - **ClinGen Allele Registry**: Canonical allele identifiers
        - **MyVariant.info**: Comprehensive variant annotations
        - **Ensembl VEP**: Variant effect predictions
        """)
        
        st.markdown("### Supported Formats")
        st.code("HGVS: NM_002496.3:c.64C>T")
        st.code("RSID: rs369602258")
        
        st.markdown("### Example Variants")
        if st.button("Load Example 1: NDUFS8"):
            st.session_state.example_input = "NM_002496.3:c.64C>T"
        if st.button("Load Example 2: BRCA1"):
            st.session_state.example_input = "NM_007294.3:c.5266dupC"
    
    # Main input
    st.markdown('<div class="section-header">Variant Input</div>', unsafe_allow_html=True)
    
    # Get input value (from example or user input)
    default_value = getattr(st.session_state, 'example_input', "")
    user_input = st.text_input(
        "Enter a genetic variant (HGVS notation or RSID):",
        value=default_value,
        placeholder="e.g., NM_002496.3:c.64C>T or rs369602258"
    )
    
    # Clear the example input after it's been used
    if hasattr(st.session_state, 'example_input'):
        delattr(st.session_state, 'example_input')
    
    analyze_button = st.button("🔬 Analyze Variant", type="primary")
    
    if analyze_button and user_input:
        # Initialize router
        router = GenomicQueryRouter()
        classification = router.classify_query(user_input)
        
        if not classification.is_genomic:
            st.error("Invalid input format. Please provide a valid HGVS notation or RSID.")
            st.stop()
        
        # Display classification info
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.write(f"**Detected identifier:** {classification.extracted_identifier}")
        st.write(f"**Type:** {classification.query_type}")
        st.markdown('</div>', unsafe_allow_html=True)
        
        try:
            start_time = time.time()
            
            # Step 1: Query ClinGen
            st.markdown('<div class="section-header">Step 1: ClinGen Allele Registry</div>', unsafe_allow_html=True)
            clingen_raw = query_clingen_allele(classification.extracted_identifier)
            clingen_data = parse_caid_minimal(clingen_raw)
            
            # Display ClinGen results
            col1, col2 = st.columns(2)
            with col1:
                st.metric("CAid", clingen_data.get('CAid', 'N/A'))
                st.metric("RSID", clingen_data.get('rsid', 'N/A'))
            with col2:
                st.metric("MANE Ensembl", clingen_data.get('mane_ensembl', 'N/A')[:20] + "..." if clingen_data.get('mane_ensembl') else 'N/A')
                st.metric("MyVariant ID", clingen_data.get('myvariant_hg38', 'N/A')[:20] + "..." if clingen_data.get('myvariant_hg38') else 'N/A')
            
            # Step 2: Get annotations
            st.markdown('<div class="section-header">Step 2: Variant Annotations</div>', unsafe_allow_html=True)
            annotations = get_variant_annotations(clingen_data)
            
            # Display any errors
            if annotations['errors']:
                for error in annotations['errors']:
                    st.warning(error)
            
            # Step 3: Display results
            if annotations['myvariant_data'] or annotations['vep_data']:
                st.markdown('<div class="section-header">Step 3: Analysis Results</div>', unsafe_allow_html=True)
                
                # Create tabs for different views
                tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary", "🧬 VEP Analysis", "🏥 MyVariant/ClinVar", "📋 Raw Data"])
                
                with tab1:
                    # Basic variant info
                    st.subheader("Variant Information")
                    if annotations['myvariant_data']:
                        myv_data = annotations['myvariant_data']
                        
                        # Basic info - handle nested data structure
                        info_cols = st.columns(3)
                        
                        # Extract chromosome info
                        chrom = (myv_data.get('hg38', {}).get('chr') or 
                                myv_data.get('chr') or 
                                myv_data.get('chrom') or 'N/A')
                        
                        # Extract position info from various possible locations
                        hg38_data = myv_data.get('hg38', {})
                        pos = (hg38_data.get('start') or hg38_data.get('end') or hg38_data.get('pos') or
                              myv_data.get('pos') or myv_data.get('start') or 
                              myv_data.get('vcf', {}).get('position') or 'N/A')
                        
                        # Extract ref/alt
                        ref = (myv_data.get('hg38', {}).get('ref') or 
                              myv_data.get('ref') or 
                              myv_data.get('vcf', {}).get('ref') or 'N/A')
                        alt = (myv_data.get('hg38', {}).get('alt') or 
                              myv_data.get('alt') or 
                              myv_data.get('vcf', {}).get('alt') or 'N/A')
                        
                        with info_cols[0]:
                            st.metric("Chromosome", chrom)
                        with info_cols[1]:
                            st.metric("Position", pos)
                        with info_cols[2]:
                            st.metric("Change", f"{ref} → {alt}")
                        
                        # Gene info
                        gene_name = (myv_data.get('genename') or 
                                   myv_data.get('gene') or 
                                   myv_data.get('symbol') or 'N/A')
                        if gene_name != 'N/A':
                            st.metric("Gene", gene_name)
                
                with tab2:
                    # Clinical significance
                    st.subheader("Clinical Significance")
                    if annotations['myvariant_data']:
                        clinical_info = display_clinical_significance(annotations['myvariant_data'])
                        if clinical_info:
                            for key, value in clinical_info.items():
                                st.write(f"**{key}:** {value}")
                        else:
                            st.info("No ClinVar clinical significance data available")
                    
                    # Functional predictions
                    st.subheader("Functional Predictions")
                    if annotations['myvariant_data']:
                        predictions = display_prediction_scores(annotations['myvariant_data'])
                        if predictions:
                            pred_cols = st.columns(2)
                            col_index = 0
                            for key, value in predictions.items():
                                with pred_cols[col_index % 2]:
                                    # Handle different value types safely
                                    try:
                                        if isinstance(value, (int, float)):
                                            st.metric(key, f"{value:.3f}" if isinstance(value, float) else str(value))
                                        else:
                                            st.metric(key, str(value))
                                    except Exception as e:
                                        st.write(f"**{key}:** {value}")
                                col_index += 1
                        else:
                            st.info("No functional prediction scores available")
                
                with tab3:
                    # Population frequencies
                    st.subheader("Population Frequencies")
                    if annotations['myvariant_data']:
                        freq_data = create_frequency_chart(annotations['myvariant_data'])
                        if freq_data is not None:
                            st.bar_chart(freq_data.set_index('Population')['Frequency'])
                        else:
                            st.info("No population frequency data available")
                
                with tab4:
                    # Raw data
                    st.subheader("Raw API Responses")
                    
                    # ClinGen data
                    with st.expander("ClinGen Allele Registry Data"):
                        st.json(clingen_data)
                    
                    # MyVariant data
                    if annotations['myvariant_data']:
                        with st.expander("MyVariant.info Data"):
                            st.json(annotations['myvariant_data'])
                    
                    # VEP data
                    if annotations['vep_data']:
                        with st.expander("Ensembl VEP Data"):
                            st.json(annotations['vep_data'])
                    
                    # Download buttons
                    st.subheader("Download Data")
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        if st.download_button(
                            "Download ClinGen Data",
                            json.dumps(clingen_data, indent=2),
                            file_name=f"clingen_{classification.extracted_identifier}.json",
                            mime="application/json"
                        ):
                            st.success("ClinGen data downloaded!")
                    
                    with col2:
                        if annotations['myvariant_data'] and st.download_button(
                            "Download MyVariant Data",
                            json.dumps(annotations['myvariant_data'], indent=2),
                            file_name=f"myvariant_{classification.extracted_identifier}.json",
                            mime="application/json"
                        ):
                            st.success("MyVariant data downloaded!")
                    
                    with col3:
                        if annotations['vep_data'] and st.download_button(
                            "Download VEP Data",
                            json.dumps(annotations['vep_data'], indent=2),
                            file_name=f"vep_{classification.extracted_identifier}.json",
                            mime="application/json"
                        ):
                            st.success("VEP data downloaded!")
            
            # Processing time
            processing_time = time.time() - start_time
            st.success(f"Analysis completed in {processing_time:.2f} seconds")
            
        except Exception as e:
            st.error(f"Analysis failed: {str(e)}")
            st.exception(e)
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style="text-align: center; color: #666;">
        <p>Built with Streamlit • Data from ClinGen, MyVariant.info, and Ensembl VEP</p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
