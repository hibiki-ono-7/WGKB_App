import streamlit as st
import pandas as pd
from pycirclize import Circos
import matplotlib.pyplot as plt
from ncbi.datasets import GeneApi, GenomeApi
from ncbi.datasets.openapi import ApiClient
from ncbi.datasets.openapi.rest import ApiException
import numpy as np
from Bio import Entrez
import textwrap
import time
        
def main() -> None:

    st.title('Custom Circos Plot Generator with Gene Metadata Integration')
    st.markdown("###### Created by Adam Hetherwick, Hibiki Ono, Nitin Kanchi, and Dr. Paulo Zaini")

    # Feedback button
    st.sidebar.markdown("### Feedback & Suggestions")
    if st.sidebar.button("Report Feedback"):
        feedback_url = "https://docs.google.com/forms/d/e/1FAIpQLSc_PnUobRzr_UXhUo_tJB2HTe4m-qJONoN1K1oPMFgTatUsqA/viewform?usp=header"
        st.markdown(f"[Click here to report feedback or issues!]({feedback_url})")
        st.stop()

    # Allow user to select which species they would like to visualize
    species_selection = st.selectbox('Select the genome you would like to visualize', ['Juglans regia', 'Juglans microcarpa x Juglans regia', 'Other'])
    
    if species_selection == 'Other':
        species_query = st.text_input('Enter the name of the species you would like to visualize.')
        
        if species_query:
            suggestions = fetch_species_suggestions(species_query)
            if suggestions:
                species_selection = st.selectbox('Did you mean:', suggestions)
            else:
                st.warning('No suggestions found. Please check the spelling.')
        st.markdown('To find the gene metadata file:\n'
                    '1. Go to https://www.ncbi.nlm.nih.gov/ \n'
                    '2. Search your species (Prunus persica for example) \n'
                    '3. Select your species then scroll down and click \"View annotated genes\" \n'
                    '4. On the table, click the square next to \'Genomic Location\' to select all rows \n'
                    '5. Click download as table, one sequence per gene')
        url = st.file_uploader('Upload genome annotation file with gene list (must be .txt or .tsv)', type=["txt", "tsv"])

    if species_selection:
        data = get_chrom_locations(species_selection)
    
    if species_selection == 'Juglans regia':
        url = 'https://raw.githubusercontent.com/nitink23/WGKB/main/juglans_regia.tsv'

    elif species_selection == 'Juglans microcarpa x Juglans regia':
        url = 'https://raw.githubusercontent.com/nitink23/WGKB/main/juglans_microcarpa.tsv'

    st.markdown("#### Enter Gene IDs and/or a gene expression file to visualize")

    # Allow user to upload optional gene expression file
    gene_exp_file = st.file_uploader('Upload a gene expression file (optional, must be .csv, .xls or .xlsx)', type=["csv", "xls", "xlsx"])

    # Read genome metadata file straight from GitHub or user uploaded file
    if url:
        genome_meta = pd.read_csv(url, delimiter='\t')

    full_data, full_data_cols = None, []

    if gene_exp_file:
        try:
            gene_exp_df = read_gene_exp_file(gene_exp_file)
            gene_id_col = st.selectbox('Select which column has the Gene IDs', gene_exp_df.columns)
            full_data = pd.merge(gene_exp_df, genome_meta, left_on=gene_id_col, right_on='Gene ID', how='inner')
            full_data_cols = full_data.columns

        except KeyError:
            st.warning('WARNING: The columns entered cannot be parsed correctly.')

    # Initialize session state for past searches
    st.session_state.setdefault("past_gene_ids", [])

    # Dropdown to show past searches
    if st.session_state.past_gene_ids:
        st.markdown(f"Gene IDs from past searches: {[int(gene_id.strip()) for gene_id in st.session_state.past_gene_ids if gene_id]}")

    # User input for multiple gene IDs
    gene_id_input = st.text_input('Enter Gene IDs (space-separated) to visualize individually.')
    genomic_ranges = []

    if gene_id_input:
            
        gene_ids = [gene_id.strip() for gene_id in gene_id_input.split() if gene_id.strip().isdigit()]

        # Save unique inputs to session state
        for gene_id in gene_ids:
            if gene_id not in st.session_state.past_gene_ids:
                st.session_state.past_gene_ids.append(gene_id)
        
        # Ensure gene_ids are in correct format
        for gene_id in gene_ids:
            for char in gene_id:
                if not char.isdigit():
                    st.error(f'ERROR: {gene_id} is not a valid Gene ID')
                    return

        gene_metadata = fetch_gene_metadata(gene_ids)
        genomic_ranges = display_gene_meta(gene_metadata, gene_ids)

        # Allow user to select a color for each gene ID
        for index, row in enumerate(genomic_ranges):

            # Assign a unique key to each color picker using the gene ID and index
            color = st.color_picker(f"Pick a color for Gene ID {row[4]}", '#ff0000', key=f"color_picker_{row[4]}_{index}")
            row.append(color)

        # Check to see if genomic ranges are valid, if invalid display warning and stop function
        if len(genomic_ranges) <= 0:
            gene_id_input = []
            genomic_ranges = []
            st.error(f"No valid genomic ranges found for Gene IDs: {', '.join(gene_ids)}")
            return

        if gene_exp_file:

            # Filter data
            selected_data = full_data[full_data['Gene ID'].astype(str).str.strip().isin(gene_ids)]


            if selected_data.empty:
                st.warning(
                    "No matching Gene IDs found. "
                    "Please ensure the Gene IDs match the dataset. "
                    "Check for data type mismatches or extra spaces."
                )
            else:
                st.write("Matching Data:", selected_data)

    show_example_plots_bool = True
    z_score_rem = 4

    if full_data is not None or gene_id_input:
        
        # Allow user to specify how many tracks they want to visualize
        available_tracks = [1, 2, 3, 4, 5] if full_data is not None else [1]
        num_tracks = st.selectbox("How many tracks would you like to visualize? (Upload gene expression file to view more than 1 track)", available_tracks)
        track_correction = 0
        line, bar = False, False

        if full_data is not None and num_tracks >= 1:
            st.markdown('Bar and line tracks can only be visualized on bottom track.')
            bar = st.checkbox("Would you like one of these tracks to be a bar plot?")
            line = st.checkbox("Would you like one of these tracks to be a line plot?")
            track_correction +=1 if bar else 0
            track_correction += 1 if line else 0
        
        # For each track, allow user to specify what data goes on which track
        track_cols = []
        bar_color = None
        line_color = None

        # Add 'Gene Location' column if user entered Gene IDs
        available_cols = list(full_data_cols)
        if genomic_ranges:
            available_cols = ['Gene Location'] + available_cols

        # Allow user to select which data they want in which tracks
        for track in range(num_tracks-track_correction):

            desired_col = st.selectbox(f"Select the data you would like to visualize in track {track+1}:", available_cols)
            include_gene_loc = False

            if desired_col == 'Gene Location':
                
                desired_data = genomic_ranges
                omit_pct = 0
                rem_outliers_bool = False
                z_score_rem = 4

            else:

                if gene_id_input:
                    include_gene_loc = st.checkbox('Would you like to highlight selected genes in this track?', key='include_gene_loc_' + str(track))

                desired_data = full_data[desired_col]

                rem_outliers_bool = st.checkbox('Would you like to remove outliers for this track?', key=f'outlier_removal_{track}')

                if rem_outliers_bool:
                    z_score_rem = st.slider(label='Choose a z-score cutoff to exclude outliers',
                                            min_value=0.1,
                                            max_value=3.5,
                                            value=3.5,
                                            step=0.01,
                                            key='track_outlier_slider_' + str(track))
                

                if invalid_col(full_data, desired_col):
                    return
                
                if st.checkbox('Would you like to omit points close to the median?', key=f'omit_points_{track}'):

                    # Add slider for user to omit certain values close to median
                    omit_pct = st.slider(
                        label="Choose a cutoff to omit points close to the median",
                        min_value=0,
                        max_value=100,
                        value=0,
                        step=1,
                        key='track_slider_' + str(track)
                    )
                else:
                    omit_pct = 0

            track_cols.append([desired_col, 'dot', desired_data, omit_pct, rem_outliers_bool, z_score_rem, include_gene_loc])
        
        # Add data to track_cols dict when user wants line or bar track
        for track_type in ['line' if line else None, 'bar' if bar else None]:

            if track_type:

                # Line and bar tracks should not visualize Gene Location
                if 'Gene Location' in available_cols:
                    available_cols.remove('Gene Location')

                desired_col = st.selectbox(f"Select which data you would like to visualize in {track_type} track:", available_cols)
                include_gene_loc = False

                if track_type == 'bar':
                    bar_color = st.color_picker(f"Pick a color for bar track", '#ff0000', key="color_picker_bar_track")

                elif track_type == 'line':
                    line_color = st.color_picker(f'Pick a color for line track', '#ff0000', key='color_picker_line_track')

                if invalid_col(full_data, desired_col):
                    return
                
                if desired_col != 'Gene Location':
                    
                    desired_data = full_data[desired_col]

                    if gene_id_input:
                        include_gene_loc = st.checkbox('Would you like to highlight selected genes in this track?', key='include_gene_loc_' + str(track_type))

                    rem_outliers_bool = st.checkbox('Would you like to remove outliers for this track?', key='outlier_removal_' + str(track_type))

                    if rem_outliers_bool:
                        z_score_rem = st.slider(label='Choose a z-score cutoff to exclude outliers',
                                                min_value=0.1,
                                                max_value=3.5,
                                                value=3.5,
                                                step=0.01,
                                                key='track_outlier_slider_' + str(track_type))

                    track_cols.append([desired_col, track_type, desired_data, 0, rem_outliers_bool, z_score_rem, include_gene_loc])

        # Plot genes not associated with any chromosome
        if full_data['Chromosome'].isna().any():

            st.markdown(f"\nThe annotated genes table for {species_selection} has {full_data['Chromosome'].isna().sum()} genes not associated with a chromosome.")

            if st.checkbox('Would you like to plot the genes not associated with a chromosome? \n'
                           'This will be a separate plot where the x-axis is index in the table, and y-axis is a user selected column.'):

                desired_col = st.selectbox(f"Select which data you would like to visualize in remaining genes track:", available_cols)
                rem_genes_color = st.color_picker(f"Pick a color for remaining genes plot", '#ff0000', key=f"rem_genes_color")
                omit_pct = st.slider(
                        label="Choose a cutoff to omit points close to the median",
                        min_value=0,
                        max_value=100,
                        value=0,
                        step=1,
                        key='track_slider_rem_genes'
                        )
                
                rem_outliers_bool = st.checkbox('Would you like to remove outliers for this track?', key='outlier_removal_rem_genes_track')

                if rem_outliers_bool:
                    z_score_rem = st.slider(label='Choose a z-score cutoff to exclude outliers',
                                            min_value=0.1,
                                            max_value=3.5,
                                            value=3.5,
                                            step=0.01,
                                            key='track_outlier_slider_rem_genes_track')

                if st.button('Click to plot genes without chromosome'):
                    if desired_col != 'Gene Location':
                        plot_rem_genes(full_data, desired_col, omit_pct, rem_genes_color, species_selection, rem_outliers_bool, z_score_rem)
                        show_example_plots_bool = False
                    else:
                        st.warning('Cannot visualize \'Gene Location\' column on remaining genes track.')

        # Plot main Circos plot
        try:
            if st.button('Click to plot!'):

                st.write('Plotting!')
                display_circos_plot(data, full_data, track_cols, bar_color, line_color, genomic_ranges, species_selection, genome_meta)
            
                show_example_plots_bool = False

        except (KeyError):
            st.error('WARNING: There was an error displaying the plot.')
            return
    
    if show_example_plots_bool:
        show_example_plots()


def fetch_species_suggestions(query):
    """Fetch species suggestions using the Entrez API."""
    try:
        handle = Entrez.esearch(db="taxonomy", term=query, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.4)  # Wait 400 ms (3 requests/second = 1 request/0.333s)
        
        if 'IdList' in record and record['IdList']:
            species_names = []
            for tax_id in record['IdList']:
                tax_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                tax_record = Entrez.read(tax_handle)
                tax_handle.close()
                species_names.append(tax_record[0]['ScientificName'])
                time.sleep(0.4)  # Rate limiting for subsequent requests
            return species_names
        else:
            return []
    except Exception as e:
        st.error(f"Error fetching species suggestions: {e}")
        return []


def show_example_plots() -> None:
    st.markdown("#### Example Plots")
    example_plot_urls = [
        "https://raw.githubusercontent.com/AJHetherwick/WGKB_App/9b5bdbadcae8e9ebcc9d4a284820d0ed19fbe9fb/juglans_regia_plot.png",
        "https://raw.githubusercontent.com/AJHetherwick/WGKB_App/9b5bdbadcae8e9ebcc9d4a284820d0ed19fbe9fb/juglans_microcarpa_plot.png",
        "https://raw.githubusercontent.com/AJHetherwick/WGKB_App/9b5bdbadcae8e9ebcc9d4a284820d0ed19fbe9fb/juglans_microcarpa_rem_genes.png"
    ]

    for url in example_plot_urls:
        st.image(url, use_column_width=True)


def get_chrom_locations(organism_name: str) -> pd.DataFrame:

    client = ApiClient()
    genome_api = GenomeApi(client)

    # Search for genome assemblies by organism name
    genome_summaries = genome_api.assembly_descriptors_by_taxon(organism_name)
    chromosomes = []
    sizes = []

    if genome_summaries.assemblies:
        for chrom_dict in genome_summaries.assemblies[0]['assembly']['chromosomes']:
            if chrom_dict.accession_version != None:
                chromosomes.append(chrom_dict.accession_version)
                sizes.append(int(chrom_dict.length))
            
    elif (not chromosomes) or (not sizes):
        st.warning(f"No genome found for organism: {organism_name}")

    return pd.DataFrame({'Chromosome': chromosomes, 'Size (bp)': sizes})


def fetch_gene_metadata(gene_ids):
    # Fetch gene metadata for all gene IDs
    api_client = ApiClient()  # Properly initialize the API client
    gene_api = GeneApi(api_client)
    try:
        gene_metadata = gene_api.gene_metadata_by_id([int(gene_id) for gene_id in gene_ids])
        return gene_metadata
    except ApiException as e:
        st.error(f"An error occurred while fetching gene metadata: {e}")
        return None


def invalid_col(full_data, desired_col) -> bool:
    # Checks for string values in given column
    if desired_col != 'Gene Location':
        if full_data[desired_col].apply(lambda x: isinstance(x, str)).all():
            st.error(f'The \'{desired_col}\' column contains string values and thus cannot be visualized on the Circos plot. The app can only read integers or floats.')
            return True
    return False


def display_formatted_response(response):
    # Function to format and display the full API response as a useful table
    if response and 'genes' in response:
        gene_info_list = []

        # Iterate through each gene in the response
        for gene_data in response['genes']:
            gene_info = gene_data.get('gene', {})
            annotations = gene_info.get('annotations', [])
            genomic_ranges = gene_info.get('genomic_ranges', [])
            chromosomes = ', '.join(gene_info.get('chromosomes', []))
            symbol = gene_info.get('symbol', 'N/A')
            gene_id = gene_info.get('gene_id', 'N/A')
            description = gene_info.get('description', 'N/A')
            taxname = gene_info.get('taxname', 'N/A')
            common_name = gene_info.get('common_name', 'N/A')

            for annotation in annotations:
                for genomic_range in genomic_ranges:
                    for loc in genomic_range.get('range', []):
                        accession_version = genomic_range.get('accession_version', 'N/A')
                        begin = loc.get('begin', 'N/A')
                        end = loc.get('end', 'N/A')
                        orientation = loc.get('orientation', 'N/A')

                        # Append gene metadata to the list
                        gene_info_list.append({
                            'Gene ID': gene_id,
                            'Symbol': symbol,
                            'Description': description,
                            'Taxname': taxname,
                            'Common Name': common_name,
                            'Chromosomes': chromosomes,
                            'Accession Version': accession_version,
                            'Range Start': begin,
                            'Range End': end,
                            'Orientation': orientation,
                            'Release Date': annotation.get('release_date', 'N/A'),
                            'Release Name': annotation.get('release_name', 'N/A')
                        })

        # Convert the list of gene info to a DataFrame
        gene_df = pd.DataFrame(gene_info_list)

        # Adjust the index to start from 1 instead of the default 0
        gene_df.index = range(1, len(gene_df) + 1)  # This sets the index starting from 1

        # Display the DataFrame with the adjusted index
        st.dataframe(gene_df)

    else:
        st.warning("No gene information available in the API response.")


def read_gene_exp_file(gene_exp_file):
    
    # See if column column headers are in row 1 or row 2
    col_row = st.selectbox('Select which row the column headers are in:', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    # Re-read gene_exp_df with correct headers
    if str(gene_exp_file.name).endswith('.csv'):
        gene_exp_df = pd.read_csv(gene_exp_file, header = col_row-1)
    else:
        gene_exp_df = pd.read_excel(gene_exp_file, header = col_row-1)

    return gene_exp_df
    

def get_colormap(data_col):
    # Convert data_col to a numpy array
    vectorized_color_col = np.array(data_col)

    # Handle case where data_col is empty or has NaNs
    if vectorized_color_col.size == 0 or np.isnan(vectorized_color_col).all():
        return ['#ff0000'] * len(data_col)
    
    # Filter out NaNs for normalization
    finite_data = vectorized_color_col[np.isfinite(vectorized_color_col)]
    min_val = finite_data.min() if finite_data.size > 0 else 0
    max_val = finite_data.max() if finite_data.size > 0 else 1
    
    # Normalize the array, ensuring the range is 0-1
    if min_val == max_val:
        normalized_arr = np.zeros_like(vectorized_color_col)  # Set all to 0 if min and max are equal
    else:
        normalized_arr = (vectorized_color_col - min_val) / (max_val - min_val)

    # Apply colormap
    cmap = plt.get_cmap('coolwarm')
    colors = [cmap(value) if np.isfinite(value) else (0.5, 0.5, 0.5, 1) for value in normalized_arr]

    return colors


def display_gene_meta(gene_metadata, gene_ids: list) -> list:
    # Log the entire API response as a table
    st.write("Full API Response:")
    display_formatted_response(gene_metadata.to_dict())  # Convert to dict for tabular format

    # Check if genes exist in the response
    if hasattr(gene_metadata, 'genes') and len(gene_metadata.genes) > 0:
        st.write(f"Fetched metadata for Gene IDs: {', '.join(gene_ids)}")

        # Parse the genomic ranges including orientation
        genomic_ranges = []
        for gene_data in gene_metadata.genes:
            gene_info = gene_data.gene

            if 'genomic_ranges' in gene_info and gene_info['genomic_ranges']:
                for genomic_range in gene_info['genomic_ranges']:
                    chrom = genomic_range['accession_version']
                    for loc in genomic_range['range']:
                        start = int(loc['begin'])
                        end = int(loc['end'])
                        orientation = loc.get('orientation', 'plus')  # Get 'orientation' from the correct location
                        genomic_ranges.append([chrom, start, end, orientation, gene_info['gene_id']])
            else:
                st.warning(f"No genomic ranges available for gene {gene_info['gene_id']}")
        return genomic_ranges
    else:
        st.warning(f'No gene ID metadata found for {gene_ids}')


def get_chrom_num(key: str, genome_meta=None) -> str:

    chrom_name = None

    # If no match in pre-defined dictionaries, check the metadata
    if genome_meta is not None:
        if 'Accession' in genome_meta.columns and 'Chromosome' in genome_meta.columns:
            chrom_dict_others = dict(zip(
                genome_meta['Accession'],
                genome_meta['Chromosome']
            ))
            chrom_name = chrom_dict_others.get(key)

    # Fallback to key if no match found
    return chrom_name or key


def omit_median_prox_data(data, omit_pct: int):

    lower_bound = data.quantile(omit_pct / 200)
    upper_bound = data.quantile(1 - (omit_pct / 200))

    return (data < lower_bound) | (data > upper_bound)


def display_circos_plot(data: dict, full_data, track_cols: list, bar_color, line_color, genomic_ranges, species_selection, genome_meta) -> None:
    
    # Get desired_data maximum to give space for yticks
    maxes = [max(desired_data) for desired_col, _, desired_data, _, _, _, _ in track_cols if desired_col != 'Gene Location']

    # Prepare Circos plot with sectors (using chromosome sizes)
    sectors = {str(row[1]['Chromosome']): row[1]['Size (bp)'] for row in data.iterrows()}
    circos = Circos(sectors, space=4, start=0, end=360 - (len(str(round(max(maxes))))+5 if maxes else 0))

    lower = 105

    # Get full data subset of the genes selected
    genomic_ranges_dict = {int(id): color for _, _, _, _, id, color in genomic_ranges}

    # Filter `full_data` by `genomic_ranges_ids`
    genomic_ranges_ids = list(genomic_ranges_dict.keys())

    # Format: track_cols[0] = [desired_col, plot_type, desired_data, omit_pct, include_gene_loc]
    for index, [desired_col, plot_type, desired_data, omit_pct, rem_outliers_bool, z_score_rem, include_gene_loc] in enumerate(track_cols):

        # omit_pct correction
        omit_pct = np.abs(100-omit_pct)

        upper = lower - 5
        if desired_col == 'Gene Location':
            width = 5
        elif plot_type == 'bar':
            width = 75/len(track_cols)
        else:
            width = 50/len(track_cols)
        lower = upper - width

        # Remove outliers for each chromosome
        if rem_outliers_bool:

            def remove_outliers(group):
                z_scores = (group[desired_col] - group[desired_col].mean()) / group[desired_col].std()
                return group[abs(z_scores) < z_score_rem]

            sector_df = (
                full_data.groupby('Chromosome', group_keys=False)
                .apply(remove_outliers)
                .reset_index(drop=True)
            )

        else:
            sector_df = full_data

        gene_full_data = sector_df[sector_df['Gene ID'].isin(genomic_ranges_ids)].copy()
        gene_full_data['gene_color'] = gene_full_data['Gene ID'].map(genomic_ranges_dict)

        # Check to see if user wants to visualize genes that are excluded because they are outliers
        for gene_id in genomic_ranges_ids:
            if int(gene_id) not in gene_full_data['Gene ID'].astype(int).values:
                st.warning(f'Gene ID {gene_id} has a value considered to be an outlier and will be excluded from {desired_col} track.')

        for chrom_num, sector_obj in enumerate(circos.sectors):
            # Map the sector name using get_chrom_num
            sector_name = get_chrom_num(sector_obj.name, genome_meta)

            if sector_name is None:
                st.warning(f"Warning: No mapping found for sector: {sector_obj.name}")
                
            # Plot sector name
            sector_obj.text(f"{sector_name}", r=110, size=10)
    
            # Add sector with correct sizing per number of tracks
            track = sector_obj.add_track((lower, upper), r_pad_ratio=0.1)
            track.axis()

            # Add x-ticks for gene location
            if index == 0:  # Apply ticks only to the outermost track

                # Major ticks (every 20 Mb)
                track.xticks_by_interval(
                    interval=20 * 1_000_000,
                    label_orientation="vertical",
                    tick_length=2,
                    label_formatter=lambda v: f"{v / 1_000_000:.0f}Mb",
                    label_size=5
                )

                # Minor ticks (every 5 Mb)
                track.xticks_by_interval(
                    interval=5 * 1_000_000,
                    tick_length=1,
                    show_label=False,  # No labels for minor ticks
                )
                
            if desired_col == 'Gene Location':
                
                # Add y-ticks only on the left side of the chr01 sector
                if chrom_num == 0:  
                    track.yticks([0, 1], list("-+"), vmin=0, vmax=1, side="left")

            # Only add tick marks to tracks that are not bar or line    
            elif plot_type not in ['bar', 'line']:
                if chrom_num == 0: 
                    num_ticks = 5 if len(track_cols) <= 3 else 3
                    track.yticks(y=np.linspace(sector_df[desired_col].min(), sector_df[desired_col].max(), num=num_ticks), \
                                labels=[f"{round(tick, 1)}" for tick in np.linspace(sector_df[desired_col].min(), sector_df[desired_col].max(), num=num_ticks)], \
                                vmin=sector_df[desired_col].min(), vmax=sector_df[desired_col].max(), side="left", label_size=7-index)
        
            # If given track is not Gene Location
            if desired_col != 'Gene Location':

                # Get subset from full_data of all the rows with a specific chromosome
                chr_data = sector_df[sector_df['Chromosome'] == sector_name].reset_index(drop=True)

                if full_data is not None:

                    gene_chr_data = gene_full_data[gene_full_data['Chromosome'] == str(sector_name)].reset_index(drop=True)

                    if plot_type == 'dot':

                        # Remove points close to the median if user selected
                        if omit_pct > 0:

                            chr_data = chr_data[omit_median_prox_data(data=chr_data[desired_col], omit_pct=omit_pct)]

                        # Get colors for column
                        chr_data['color_' + str(index)] = get_colormap(chr_data[desired_col])

                        track.scatter(chr_data['Begin'].tolist(), chr_data[desired_col].tolist(), color=chr_data['color_' + str(index)], 
                                      vmin=sector_df[desired_col].min(), vmax=sector_df[desired_col].max(), s=5)

                        # If user wants selected gene to be highlighted on current track
                        if include_gene_loc:

                            track.scatter(gene_chr_data['Begin'].tolist(), gene_chr_data[desired_col].tolist(), color=gene_chr_data['gene_color'], 
                                          vmin=sector_df[desired_col].min(), vmax=sector_df[desired_col].max(), s=15)

                    if not pd.isna(chr_data[desired_col].max()):

                        if plot_type == 'line':

                            # Add solid line connecting data points
                            chr_data_sorted = chr_data.sort_values(by='Begin')
                            track.line(chr_data_sorted['Begin'].tolist(), chr_data_sorted[desired_col].tolist(), color=line_color, 
                                       vmin=sector_df[desired_col].min(), vmax=sector_df[desired_col].max())

                            # If user wants selected gene to be highlighted on current track
                            if include_gene_loc:

                                track.scatter(gene_chr_data['Begin'].tolist(), gene_chr_data[desired_col].tolist(), color=gene_chr_data['gene_color'], 
                                              vmin=sector_df[desired_col].min(), vmax=sector_df[desired_col].max(), s=15, zorder=2)

                        elif plot_type == 'bar':

                            # Ensure the column has values in it
                            vmin = 0 if sector_df[desired_col].min() > 0 else sector_df[desired_col].min()
                            track.bar(chr_data['Begin'].tolist(), chr_data[desired_col].tolist(), ec=bar_color, lw=0.9, vmin=vmin, vmax=sector_df[desired_col].max())

                            # If user wants selected gene to be highlighted on current track
                            if include_gene_loc:
                                track.scatter(gene_chr_data['Begin'].tolist(), gene_chr_data[desired_col].tolist(), color=gene_chr_data['gene_color'], 
                                             vmin=vmin, vmax=sector_df[desired_col].max(), s=15)

            # If user manually included gene IDs
            else:
                for chrom, x, _, orientation, gene_id, color_to_use in desired_data:
                    if chrom == sector_obj.name:
                        y = 1 if orientation.value == 'plus' else 0
                        track.scatter([x], [y], color=color_to_use, s=20)

    # Add colorbar for expression level columns
    exp_cols = [[desired_col, plot_type, index] for index, [desired_col, plot_type, _, _, _, _, _] in enumerate(track_cols) if desired_col != 'Gene Location']

    for index, (label, plot_type, track_index) in enumerate(exp_cols):

        if plot_type == 'bar' or plot_type == 'line':

            # Shorten title if we are examining Expression level
            if label == 'Expression level (average normalized FPKM of all 36 samples)':
                label = 'Average normalized FPKM expression'

            circos.colorbar(
                bounds=(0.92, 1+index*0.1, 0.2, 0),
                vmin=sector_df[label].min(),
                vmax=sector_df[label].max(),
                orientation="horizontal",
                label=f'Track {track_index+1}: ' + label,
                label_kws=dict(size=8, color="black"),
                tick_kws=dict(labelsize=8, colors="black")
            )

        else:
            circos.colorbar(
                bounds=(0.92, 1+index*0.1, 0.2, 0.02),
                vmin=sector_df[label].min(),
                vmax=sector_df[label].max(),
                cmap="coolwarm",
                orientation="horizontal",
                label=f'Track {track_index+1}: ' + label,
                label_kws=dict(size=8, color="black"),
                tick_kws=dict(labelsize=8, colors="black")
            )

    plt.tight_layout()

    # Render the plot using Matplotlib
    fig = circos.plotfig()
    
    add_species_title(species_selection, circos, num_tracks=len(track_cols))

    # Add legend for selected gene IDs
    if genomic_ranges:

        scatter_legend = circos.ax.legend(
            handles=[plt.Line2D([0], [0], color=row[-1], marker='o', ls='None', ms=8) for row in genomic_ranges],  
            labels=[row[4] for row in genomic_ranges],  
            bbox_to_anchor=(-0.12, 1.1),
            loc='upper left',
            fontsize=8,
            title="Gene ID",
            handlelength=2
        )

        scatter_legend._legend_title_box._text_pad = 5
        circos.ax.add_artist(scatter_legend)

    # Add legend for sliders 
    removal_legend_items = [
        (f"Track {track_index + 1}: Omit {omit_pct}% around median", track_index)
        for track_index, (_, _, _, omit_pct, _, _, _) in enumerate(track_cols)
        if omit_pct > 0
    ]

    # Add removal legend if there are items
    if removal_legend_items:

        # Create handles and labels for all tracks in one go
        legend_handles = [
            plt.Line2D([0], [0], color="black", linestyle='None')
            for _ in removal_legend_items
        ]
        legend_labels = [label for label, _ in removal_legend_items]

        add_median_legend(circos, legend_handles, legend_labels, x=-0.14, y=0)
        
    # Display the plot in Streamlit
    st.pyplot(fig)


def add_species_title(species_selection, circos_obj, num_tracks) -> None:
    # Add species legend in the middle of plot

    wrapped_labels = [
        "\n".join(textwrap.wrap(species_selection, width=20, break_long_words=False))
    ]

    species_legend = circos_obj.ax.legend(
        handles=[plt.Line2D([0], [0], color="black", linestyle='None') for _ in species_selection],
        labels=wrapped_labels,  # Use wrapped labels
        bbox_to_anchor=(0.48, 0.5),
        loc='center',
        prop={'size': 12 - (num_tracks)//2, 'style': 'italic'}
    )

    circos_obj.ax.add_artist(species_legend)


def add_median_legend(circos_obj, handles, labels, x: int, y: int) -> None:

    removal_legend = circos_obj.ax.legend(
        handles=handles,  
        labels=labels, 
        bbox_to_anchor=(x, y),
        title="Median Percentile Removal",
        fontsize=8,
        title_fontsize=9,
        handlelength=1.5
    )

    # removal_legend._legend_title_box._text_pad = 5
    circos_obj.ax.add_artist(removal_legend)


def plot_rem_genes(full_data, desired_col, omit_pct, color, species_selection, rem_outliers_bool, z_score_rem) -> None:

    rem_genes_df = full_data[pd.isna(full_data['Chromosome'])].reset_index()
    desired_data = rem_genes_df[desired_col]

    if rem_outliers_bool:
            
        z_scores = (desired_data - desired_data.mean()) / desired_data.std()
        sector_data = desired_data[abs(z_scores) < z_score_rem]

    else:
        sector_data = desired_data

    sectors = {rem_genes_df['Accession'].iloc[0]: max(rem_genes_df.index.tolist())}
    circos = Circos(sectors, space=4, start=0, end=360 - (len(str(round(max(sector_data))))+5))

    for sector_obj in circos.sectors:

        track = sector_obj.add_track((70, 100), r_pad_ratio=0.1)
        track.axis()

        # Add y-ticks only on the left side
        track.yticks(y=np.linspace(min(sector_data), max(sector_data), num=5), \
                     labels=[f"{round(tick, 1)}" for tick in np.linspace(min(sector_data), max(sector_data), num=5)], \
                     vmin=min(sector_data), vmax=max(sector_data), side="left", label_size=7)

        if omit_pct > 0:

            sector_data = sector_data[omit_median_prox_data(data=sector_data, omit_pct=omit_pct)]

        rem_genes_x = sector_data.index.tolist()
        sector_data = sector_data.tolist()

        track.scatter(x=rem_genes_x, y=sector_data, color=color, 
                      vmin=min(sector_data), vmax=max(sector_data), s=12)
        
        # Add y-axis gridlines manually
        track.grid()

        # Add xticks for index
        indices = list(range(len(sector_data)))

        # Major ticks
        track.xticks_by_interval(
            interval=max(indices) // 10,
            show_label=True,
            label_orientation="vertical",
            tick_length=2,
            label_size=10,
            outer=False,
        )

        # Minor ticks
        minor_interval = 1 if max(indices) // 50 == 0 else max(indices) // 50
        track.xticks_by_interval(
            interval=minor_interval,  # Minor ticks at a smaller interval
            outer=False,
            tick_length=1,  # Shorter tick marks for minor ticks
            show_label=False,  # No labels for minor ticks
        )
        
    fig_2 = circos.plotfig()
        
    # Add sector name as a legend
    scatter_legend = circos.ax.legend(
        bbox_to_anchor=(1.05, 1.05),
        loc='upper right',
        fontsize=8,
        title=sector_obj,
        handlelength=2
    )

    scatter_legend._legend_title_box._text_pad = 5
    circos.ax.add_artist(scatter_legend)

    if omit_pct > 0:

        legend_handles = [plt.Line2D([0], [0], color="black", linestyle='None')]
        legend_labels = [f"Omit {omit_pct}% around median"]

        add_median_legend(circos, legend_handles, legend_labels, x=-0.05, y=0.1)

    add_species_title(species_selection, circos, num_tracks=1)

    st.pyplot(fig_2)


main()
