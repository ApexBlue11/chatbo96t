import streamlit as st
import requests
import pandas as pd
import json
import time
import re
from typing import Dict, Any, List, Optional
from dataclasses import dataclass
import plotly.express as px
import plotly.graph_objects as go
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

def create_frequency_plot(myvariant_data):
    """Create a population frequency visualization."""
    freq_data = []
    
    # Extract gnomAD frequencies
    for key, value in myvariant_data.items():
        if 'gnomad' in key.lower() and 'af' in key.lower() and isinstance(value, (int, float)):
            if value > 0:  # Only show non-zero frequencies
                # Clean up the key for display
                clean_key = key.replace('gnomad_', '').replace('_af', '').replace('_', ' ').title()
                freq_data.append({'Population': clean_key, 'Frequency': value})
    
    if freq_data:
        df_freq = pd.DataFrame(freq_data)
        fig = px.bar(df_freq, x='Population', y='Frequency', 
                    title='Population Frequencies (gnomAD)',
                    labels={'Frequency': 'Allele Frequency'})
        fig.update_layout(xaxis_tickangle=-45)
        return fig
    return None

def display_clinical_significance(myvariant_data):
    """Display ClinVar clinical significance information."""
    clinvar = myvariant_data.get('clinvar', {})
    if not clinvar:
        return None
    
    clinical_info = {}
    
    # Extract key clinical information
    if 'clinical_significance' in clinvar:
        clinical_info['Clinical Significance'] = clinvar['clinical_significance']
    if 'review_status' in clinvar:
        clinical_info['Review Status'] = clinvar['review_status']
    if 'conditions' in clinvar:
        conditions = clinvar['conditions']
        if isinstance(conditions, list):
            clinical_info['Associated Conditions'] = ', '.join([
                cond.get('name', str(cond)) if isinstance(cond, dict) else str(cond) 
                for cond in conditions
            ])
        else:
            clinical_info['Associated Conditions'] = str(conditions)
    
    return clinical_info

def display_prediction_scores(myvariant_data):
    """Display functional prediction scores."""
    dbnsfp = myvariant_data.get('dbnsfp', {})
    if not dbnsfp:
        return None
    
    predictions = {}
    
    # SIFT
    if 'sift' in dbnsfp:
        sift_data = dbnsfp['sift']
        if isinstance(sift_data, dict):
            if 'score' in sift_data:
                predictions['SIFT Score'] = sift_data['score']
            if 'pred' in sift_data:
                predictions['SIFT Prediction'] = sift_data['pred']
    
    # PolyPhen
    if 'polyphen2' in dbnsfp:
        pp2_data = dbnsfp['polyphen2']
        if isinstance(pp2_data, dict):
            if 'hdiv' in pp2_data and isinstance(pp2_data['hdiv'], dict):
                if 'score' in pp2_data['hdiv']:
                    predictions['PolyPhen2 HDiv Score'] = pp2_data['hdiv']['score']
                if 'pred' in pp2_data['hdiv']:
                    predictions['PolyPhen2 HDiv Prediction'] = pp2_data['hdiv']['pred']
    
    # CADD
    if 'cadd' in dbnsfp:
        cadd_data = dbnsfp['cadd']
        if isinstance(cadd_data, dict):
            if 'phred' in cadd_data:
                predictions['CADD Phred Score'] = cadd_data['phred']
    
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
                tab1, tab2, tab3, tab4 = st.tabs(["📊 Summary", "🧬 Clinical", "📈 Frequencies", "📋 Raw Data"])
                
                with tab1:
                    # Basic variant info
                    st.subheader("Variant Information")
                    if annotations['myvariant_data']:
                        myv_data = annotations['myvariant_data']
                        
                        # Basic info
                        info_cols = st.columns(3)
                        with info_cols[0]:
                            st.metric("Chromosome", myv_data.get('hg38', {}).get('chr', myv_data.get('chr', 'N/A')))
                        with info_cols[1]:
                            st.metric("Position", myv_data.get('hg38', {}).get('pos', myv_data.get('pos', 'N/A')))
                        with info_cols[2]:
                            st.metric("Change", f"{myv_data.get('ref', 'N/A')} → {myv_data.get('alt', 'N/A')}")
                        
                        # Gene info
                        if myv_data.get('genename'):
                            st.metric("Gene", myv_data['genename'])
                
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
                                    st.metric(key, value)
                                col_index += 1
                        else:
                            st.info("No functional prediction scores available")
                
                with tab3:
                    # Population frequencies
                    st.subheader("Population Frequencies")
                    if annotations['myvariant_data']:
                        freq_plot = create_frequency_plot(annotations['myvariant_data'])
                        if freq_plot:
                            st.plotly_chart(freq_plot, use_container_width=True)
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
