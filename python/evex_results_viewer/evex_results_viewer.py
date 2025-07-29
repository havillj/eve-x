import sys

# Ensure config file exists or command line chosen root directory exists

# Examples: 
# python3 evex_result_viewer.py
# python3 evex_result_viewer.py -root ./
# python3 evex_result_viewer.py -root /foo
if len(sys.argv) == 1: # If root directory not specifeid in command line
    try:
        from config import RESULTS_DIR, SEQUENCES_DIR, SPECIMEN_RESULTS_DIR, VIRUS_RESULTS_DIR, DIAGRAMS_DIR
    except:
        sys.exit("Error: Config file not found")
else:
    try:
        index = sys.argv.index("-root")
        ROOT_DIR = sys.argv[index+1]
        if ROOT_DIR[-1] != '/':
             RESULTS_DIR = ROOT_DIR + "/results/"
        else:
            RESULTS_DIR = ROOT_DIR + "results/"
        SEQUENCES_DIR = "sequences/"
        DIAGRAMS_DIR = "diagrams/"
        SPECIMEN_RESULTS_DIR = RESULTS_DIR + "specimens/"
        VIRUS_RESULTS_DIR = RESULTS_DIR + "viruses/"
    except ValueError:
        sys.exit("-root was not found")
    except IndexError:
        sys.exit("Root directory not specified")

import dearpygui.dearpygui as dpg
from pathlib import Path
from pdf2image import convert_from_path
import re
from Bio import SeqIO
import os

dir_name_finder = "(\s|\w)*\/"
pattern_dir = re.compile(dir_name_finder)

matches = pattern_dir.finditer(VIRUS_RESULTS_DIR)
VIRUS_DIR = ""
for m in matches: # Get chosen name for the virus directory in the results directory
    VIRUS_DIR = m.group()

matches = pattern_dir.finditer(SPECIMEN_RESULTS_DIR)

SPECIMEN_DIR = ""
for m in matches: #Get the chosen name for the specimen directory in the results directory
    SPECIMEN_DIR = m.group()

WHITE = (255, 255, 255, 255)
BLUE = (60, 136, 238, 255)
RED = (235, 65, 59, 255)
YELLOW = (250, 179, 64, 255)
GREEN = (100, 247, 64, 255)
GRAY = (90, 90, 90, 255)
CHAR_HEIGHT = 10
CHAR_WIDTH = 8


v_sequences_path = VIRUS_RESULTS_DIR + SEQUENCES_DIR
v_diagrams_path = VIRUS_RESULTS_DIR + DIAGRAMS_DIR
png_images_path = "/tmp/png_images"

txt_regex = ".*\.txt"
txt_pattern = re.compile(txt_regex)

fasta_regex = ".*\.fasta"
fasta_pattern = re.compile(fasta_regex)

pdf_regex = ".*\.pdf"
pdf_pattern = re.compile(pdf_regex)

cluster_regex = "C([0-9]*)(\_\|\_.*)" #Seq object omits '>' character
cluster_pattern = re.compile(cluster_regex)

cluster_changes_count = 0

# Files not to show
# bam_regex = ".*\.bam"
# tsv_regex = ".*\.tsv"
# csv_regex = ".*\.csv"
# xml_regex = ".*\.xml"

exclude_file_regex = ".*\.(bam|tsv|csv|xml)"
exclude_file_pattern = re.compile(exclude_file_regex)

class Interface():

    def __init__(self, path):

        #Convert path to dictionary
        self.file_sys_dictionary = {}
        # Omit the '/' from results directory constant
        self.file_sys_dictionary[RESULTS_DIR[:-1]] = self.file_sys_to_dictionary(path)
        
        # Colors associated to each base
        self.base_color = {"A": GREEN, 
                           "C": YELLOW, 
                           "G":RED, 
                           "T": BLUE}

        #Run the display
        self.show_window()

    def file_sys_to_dictionary(self, path):
        """
            Creates a multi-layered dictionary of the file system

            Input: 
                Path p to file or directory

            Return: 
                Dictionary version of file system
        """

        if not path.is_dir(): # If the path leads to a file
            return None #Indicator that the path reached its end
        else: # Otherwise, it's a directory
            directory = {}
            for content in path.iterdir(): # Iterate through directory content
                if not exclude_file_pattern.fullmatch(content.name): # Omit .bam, .tsv, .csv, .xml files
                    directory[content.name] = self.file_sys_to_dictionary(content)
            return dict(sorted(directory.items())) #Ensures files and diretories are sorted alphabetically

    def pdf_to_png_dpg(self, name, path_str):
        """
        Converts pdf to png if not already exisiting and 
        creates a tag for the png image

        Input: 
            name of the pdf file
            path_str to the pdf file

        Return:
            width and height of the png image
        """

        # Determine whether png version of pdf already exists before creating it
        png = f"{png_images_path}/{name}0001-1.png"
        png_path = Path(png) 
        if not png_path.exists():
            convert_from_path(f"{path_str}", fmt="png", output_folder=png_images_path, output_file=name, paths_only=True)
            # Load image and create image in dearpygui
            width, height, channels, data = dpg.load_image(png)
            with dpg.texture_registry():
                dpg.add_static_texture(width, height, data, tag=png)

            return (width, height)


    def resize_image(self, sender, app_data, user_data):
        """
        Resizes the png whenever the window containing the png is resized

        Input:
            sender = field dearpygui uses to bind item to callback to the function resize_image
            app_data = field dearpygui uses
            user_data = tuple containing the alias of the window showing the image at index 0
                        and the drawlist tag that is displaying the png at index 1
        Return:
            None
        """

        window_alias = user_data[0]
        drawlist_tag = user_data[1]
        ori_width = user_data[2]
        ori_height = user_data[3]
        width, height = dpg.get_item_rect_size(window_alias) # Get the new dimensions of the resized window
        height_image = dpg.get_item_height(drawlist_tag)
        width_image = dpg.get_item_width(drawlist_tag)
        if height != height_image:
            if height_image < height:
                scale = height / ori_height
                height_image = ori_height * scale
                width_image = ori_width * scale
        if width != width_image:
            scale = width / ori_width
            height_image = ori_height * scale
            width_image = ori_width * scale
        dpg.set_item_height(drawlist_tag, height_image) # Adjust the drawlist size based on the new window size
        dpg.set_item_width(drawlist_tag, width_image)

        # Remake the new image to adjust for window resize 
        if dpg.does_alias_exist(drawlist_tag):
            dpg.delete_item(drawlist_tag, children_only=True) # Deletes the png inside of the drawlist
        dpg.draw_image(f"{png_images_path}/{window_alias}0001-1.png", (0, 0), (width_image, height_image), uv_min=(0,0), uv_max=(1, 1), parent=drawlist_tag)


    def display_file_pdf(self, sender, app_data, user_data):
        """
        Displays the file when user wishes to view the file

        Input: 
            sender = field dearpygui uses to bind item using the callback to the function display_file_pdf()
            app_data = field dearpygui uses
            user_data = tuple containing virus name at index 0 and path to file as a string at index 1
        
        Return:
            None
        """

        # Make a png directory to store all made png; delete folder once program closes 
        alias = user_data[0]
        path_to_file = user_data[1]
        png_dir = Path(png_images_path)
        if not png_dir.exists() and not png_dir.is_dir():
            try:
                png_dir.mkdir()
            except FileExistsError:
                print(f"Directory {png_dir} already exists")
            except PermissionError:
                print(f"Permission denied: Unable to create {png_dir}")
            except Exception as e:
                print(f"An error occured: {e}")
        
        # Determine whether the window has been created or not
        if not dpg.does_alias_exist(alias):
            # Allow for png image to adjust to window resize
            width_image, height_image = self.pdf_to_png_dpg(alias, path_to_file)

            with dpg.item_handler_registry(tag=f"{alias}_resize"):
                dpg.add_item_resize_handler(callback=self.resize_image, user_data=(alias, path_to_file, width_image, height_image))


            # Make the window have a maximum width and height if image is too big
            width_window = width_image
            height_window = height_image
            if width_window > 1400:
                width_window = 1400
                height_window *= 1400/width_image
            if height_window > 750:
                height_window = 750
                width_window = 750/height_image

            # if VIRUS_RESULTS_DIR not in path_to_file: # Might have an issue if DIAGRAMS_DIR == VIRUS_DIR for whatever reason
            width_image *= 0.75
            height_image *= 0.75
            with dpg.window(label=alias, tag=alias, show=True, horizontal_scrollbar=False, width=width_window, height=height_window+40, pos=(250, 250)):
                png = f"{png_images_path}/{alias}0001-1.png"
                with dpg.drawlist(tag=path_to_file, width=width_image, height=height_image, parent=alias):
                    width, height = dpg.get_item_rect_size(alias)
                    dpg.draw_image(png, (0, 0), (width, height_image), uv_min=(0, 0), uv_max = (1, 1))
            dpg.bind_item_handler_registry(alias, f"{alias}_resize") # Bind window and item handler tags
        else: # Otherwise show the window if not already
            dpg.configure_item(alias, show=True)
            dpg.focus_item(alias)


    def display_alignment(self, seq):
        """
        Helper function to display an aligned sequence

        Input:
            Sequence object

        Return:
            None
        """

        dashes = ""
        for character in seq.seq:
            if character == "-":
                dashes += "-"
            else:
                if dashes != "":
                    dpg.add_text(default_value=dashes, color=WHITE)
                    dashes = ""
                if character in self.base_color:
                    dpg.add_text(default_value=character, color=self.base_color[character])
                else:
                    dpg.add_text(default_value=character, color=WHITE)
        if dashes != "":
            dpg.add_text(default_value=dashes, color=WHITE)
    

    def copy_seq(self, sender, app_data, user_data):
        """
        Helper function to copy sequence to clipboard

        Input:
            sender = field dearpygui uses to bind item using the callback to the function cop_seq()
            app_data = field dearpygui uses
            user_data = sequence that is to be copied to clipboard
        
        Return:
            None
        """
        
        dpg.set_clipboard_text(user_data)


    def zoom_out_fasta_aligned_viruses(self, sender, app_data, user_data):
        """
        Shows a zoom out version of the fasta aligned sequences in the virus directory represented with a 
        line showing the length of the sequence amd rectangles indicating alignment (or pressence of a base
        in the fasta file)

        Input:
            sender = field dearpygui uses to bind item using the callback to teh function zoom_out_fasta__aligned_virus()
            app_data = fiedl dearpygui uses
            user_data = file name at index 0; path to file at index 1

        Return:
            None
        """

        alias = user_data[0]
        zoom_alias = alias+"_zoom"
        if not dpg.does_alias_exist(zoom_alias):
            fasta_file = user_data[1]
            seqs = list(SeqIO.parse(fasta_file, "fasta"))
            num_seqs = len(seqs)
            if num_seqs > 15:
                num_seqs = 15
            with dpg.window(label=zoom_alias, tag=zoom_alias, show=True, horizontal_scrollbar=True, width=1500, height=(46* num_seqs), pos=(250, 250)):
                cur_row = 0
                for seq in seqs:
                    with dpg.item_handler_registry(tag=f"{fasta_file}_{seq.id}_copyzoom_{cluster_changes_count}"):
                        dpg.add_item_clicked_handler(button=dpg.mvMouseButton_Left, callback=self.copy_seq, user_data=str(seq.seq))
                    self.zoom_out_simplified_seq(fasta_file, seq, cur_row)
                    dpg.bind_item_handler_registry(f"{fasta_file}_{seq.id}_zoom_{cluster_changes_count}", f"{fasta_file}_{seq.id}_copyzoom_{cluster_changes_count}")
                    cur_row += 1
        else:
            dpg.configure_item(zoom_alias, show=True)
            dpg.focus_item(zoom_alias)
        

    def display_file_fasta_aligned_viruses(self, sender, app_data, user_data):
        """
        Displays the aligned fasta files 

        Input: 
            sender = field dearpygui uses to bind item (button) using the callback to the function display_file_fasta_aligned_viruses()
            app_data = field dearpygui uses
            user_data = tuple containing name of fasta file at index 0 and path to file as a string at index 1
        
        Return:
            None
        """

        alias = user_data[0]
        if not dpg.does_alias_exist(alias): # Make window if not already
            # Parse the fasta file
            fasta_file = f"{user_data[1]}"
            seqs = list(SeqIO.parse(fasta_file, "fasta"))
            maxIDLength = max(len(seq.id) for seq in seqs)
            
            # Create the window
            num_seqs = len(seqs)
            if num_seqs > 10:
                num_seqs = 10
            with dpg.window(label=alias, tag=alias, show=True, horizontal_scrollbar=True, width=1500, height=30 + 25*num_seqs, pos=(250, 250)):
                # Create the option to view a zoom out (or simplified) version of the fasta aligned file
                with dpg.menu_bar():
                    dpg.add_menu_item(label="Zoom out", callback=self.zoom_out_fasta_aligned_viruses, user_data=(alias, fasta_file))
                
                ref_seq = seqs[0]
                # Diplay the reference sequence at the top
                with dpg.group(horizontal=True):
                    dpg.add_button(label=f"{ref_seq.id:>{maxIDLength+1}}:", callback=self.copy_seq, user_data=ref_seq.seq)
                    with dpg.tooltip(dpg.last_item()):
                        dpg.add_text("Click to copy to clipboard")
                    with dpg.item_handler_registry(tag=f"{ref_seq.id}_{fasta_file}_location"): # Make callback when sequence is clicked on
                        dpg.add_item_clicked_handler(button=dpg.mvMouseButton_Left, callback=self.copy_seq, user_data=ref_seq.seq)
                    with dpg.group(horizontal=True, horizontal_spacing=0, tag=f"{ref_seq.id}_{fasta_file}_seq"): # Make the bases shown be horizontal with each other
                        self.display_alignment(ref_seq)
                    with dpg.tooltip(parent=f"{ref_seq.id}_{fasta_file}_seq"):
                        dpg.add_text(f"Click to copy ({ref_seq.id})")
                    dpg.bind_item_handler_registry(f"{ref_seq.id}_{fasta_file}_seq", f"{ref_seq.id}_{fasta_file}_location") # Bind sequence to copy function handler

                # Display the rest of the sequences in the file
                with dpg.child_window(width=(7*(maxIDLength+len(ref_seq.seq))+61)):
                    for seq in seqs[1:]:
                        with dpg.group(horizontal=True):
                            dpg.add_button(label=f"{seq.id:>{maxIDLength}}:", callback=self.copy_seq, user_data=seq.seq)
                            with dpg.tooltip(dpg.last_item()):
                                dpg.add_text("Click to copy to clipboard")
                            with dpg.item_handler_registry(tag=f"{seq.id}_{fasta_file}_location"): # Make callback when squence is clicked on
                                dpg.add_item_clicked_handler(button=dpg.mvMouseButton_Left, callback=self.copy_seq, user_data=seq.seq)
                            with dpg.group(horizontal=True, horizontal_spacing=0, tag=f"{seq.id}_{fasta_file}_seq"): # Make the bases shown be horizontal with each other
                                self.display_alignment(seq)
                            with dpg.tooltip(parent=f"{seq.id}_{fasta_file}_seq"):
                                dpg.add_text(f"Click to copy ({seq.id})")
                            dpg.bind_item_handler_registry(f"{seq.id}_{fasta_file}_seq", f"{seq.id}_{fasta_file}_location") # Bind sequence to copy function handler

        else: # Otherwise display the window
            dpg.configure_item(alias, show=True)
            dpg.focus_item(alias)


    def display_file_fasta_not_aligned(self, sender, app_data, user_data):
        """
        Displays fasta files that aren't aligned

        Input: 
            sender = field dearpygui uses to bind item (button) using the callback to the function display_file_fasta_not_aligned()
            app_data = field dearpygui uses
            user_data = tuple containing name of fasta file at index 0 and path to file as a string at index 1
        
        Return:
            None
        """

        alias = user_data[0]
        if not dpg.does_alias_exist(alias): # Make window if not already
            # Parse the fasta file
            fasta_file = user_data[1]
            seqs = list(SeqIO.parse(fasta_file, "fasta"))

            # Make maximum window height in case fasta file has many sequences
            num_seqs = len(seqs)
            if num_seqs > 10:
                num_seqs = 10
            
            # Used to make aligning sequence ID
            maxIDLength = max(len(seq.id) for seq in seqs)

            with dpg.window(label=alias, tag=alias, show=True, horizontal_scrollbar=True, width=1080, height=30+25*num_seqs, pos=(250, 250)):
                # Display the sequences
                for seq in seqs:
                    with dpg.group(horizontal=True):
                        with dpg.item_handler_registry(tag=f"{seq.id}_{fasta_file}_location"):
                            dpg.add_item_clicked_handler(button=dpg.mvMouseButton_Left, callback=self.copy_seq, user_data=str(seq.seq))
                        dpg.add_button(label=f"{seq.id:>{maxIDLength}}: ", callback=self.copy_seq, user_data=str(seq.seq))
                        with dpg.tooltip(dpg.last_item()):
                            dpg.add_text("Click to copy to clipboard")
                        dpg.add_text(f"{seq.seq!s}", tag=f"{seq.id}_{fasta_file}_seq")
                        with dpg.tooltip(dpg.last_item()):
                            dpg.add_text(f"Click to copy {seq.id}")
                        dpg.bind_item_handler_registry(f"{seq.id}_{fasta_file}_seq", f"{seq.id}_{fasta_file}_location")
        else: # Otherwise display the window
            dpg.configure_item(alias, show=True)
            dpg.focus_item(alias)


    def exit_app(self):
        """
        Removes the png directory that stores the png files of the pdf

        Input:
            None

        Return:
            None
        """

        p = Path(png_images_path)
        if p.exists():
            os.system(f"rm -f {png_images_path}/*")
            os.system(f"rmdir {png_images_path}")

    def update_order(self, sender, app_data, user_data):
        """
        Update the order of the clsuter groups. Once the order is performed, the window redisplays itself
        with the updated cluster groups

        Input:
            sender = field dearpygui uses to bind item using the callback to the function display_file()
            app_data = field dearpygui uses
            user_data = tuple containing the alias of the cluster at index 0; cluster fasta file path at index 1

        Return:
            None 
        """

        alias = user_data[0]
        fasta_file = user_data[1]
        seqs = list(SeqIO.parse(fasta_file, "fasta"))
        ref_seq = seqs.pop(0) # Remove the reference sequence from being sorted

        # Update the changes made by the user when editing the input text fields by calling the tag assigned to each field
        for seq in seqs:
            cluster_pattern_group = cluster_pattern.findall(seq.id)[0]
            seq_name = cluster_pattern_group[1]
            changes = dpg.get_value(seq_name)
            if dpg.does_alias_exist(f"{seq_name}_popup"):
                changes = dpg.get_value(f"{seq_name}_popup")
            if not changes.isdecimal():
                return
            seq.id = f"C{changes}{seq_name}"
        
        seqs.sort(key=lambda seq: int(cluster_pattern.findall(seq.id)[0][0])) # Sort by Cluster number
        seqs.insert(0, ref_seq) # Insert reference sequence back

        SeqIO.write(seqs, fasta_file, "fasta") # Update the cluster file
        
        # Remove the previous window and display with the changes made by the user
        dpg.delete_item(alias)
        if dpg.does_alias_exist(alias+"_zoom"):
            dpg.delete_item(alias+"_zoom")
        self.display_cluster_fasta(sender=f"{alias}_button", app_data=0, user_data=(alias, fasta_file))
        

    def change_order_popup(self, sender, app_data, user_data):
        """
        Helper function to allow the user to change the cluster number when clicking on a zoomed out sequence

        Input:
            sender = field dearpygui uses to bind item using the callback to the function change_order_popup()
            app_data = field dearpygui uses
            user_data = sequence object at index 0; file name at index 1; path to file at index 2
        Return:
            None
        """
        #TODO: Add a copy sequence button
        seq = user_data[0]
        cluster_pattern_group = cluster_pattern.findall(seq.id)[0]
        seq_cluster_num = cluster_pattern_group[0]
        seq_name = cluster_pattern_group[1]

        alias = user_data[1]
        fasta_file = user_data[2]
        with dpg.popup(parent=f"{fasta_file}_{seq.id}_zoom_{cluster_changes_count}", modal=True, mousebutton=dpg.mvMouseButton_Left): #TODO: Figure out how to make it so only one click required
            with dpg.menu_bar():
                dpg.add_menu_item(label="Save changes", callback=self.update_order, user_data=(alias, fasta_file))
            with dpg.group(horizontal=True, horizontal_spacing=0):
                dpg.add_text("C")
                dpg.add_input_text(default_value=seq_cluster_num, tag=f"{seq_name}_popup", width=20, decimal=True)
                seq_id = f"{seq_name}"
                dpg.add_text(f"{seq_id:<{len(seq.id)}}")
            dpg.add_separator()
            dpg.add_button(label="Copy Sequence", callback=self.copy_seq, user_data=str(seq.seq))


    def zoom_out_simplified_seq(self, fasta_file, seq, cur_row):
        """
        Helper function to display the sequences simplified manner with a line and 
        overlapping rectangles that indicate where bases are aligned
        
        Input:
            fasta_file = used to help with internal dearpygui tag creation; unique path to each file
            seq = sequence object to display the sequence of
            cur_row = the current row to display the sequence in
        
        Return:
            None
        """

        seq_len = len(seq.seq)
        with dpg.group(tag=f"{fasta_file}_{seq.id}_zoom_{cluster_changes_count}"):
            with dpg.drawlist(height=CHAR_HEIGHT, width=seq_len, pos=(CHAR_WIDTH, CHAR_HEIGHT*cur_row)):
                dpg.draw_line((0, CHAR_HEIGHT/2), (seq_len*CHAR_WIDTH, CHAR_HEIGHT/2), thickness=0, color=GRAY) # Draw line representing the sequence
                start_rect = 0 # Left side of the rectangle
                cur_pos = 0
                while cur_pos < seq_len:
                    if seq.seq[cur_pos] != '-': # Reading a base
                        if start_rect==0: # Hasn't marked the start of the alignment yet
                            start_rect = cur_pos
                    else: #cur_pos is at a dash
                        if start_rect!=0: # Found the end of the rectangle
                            dpg.draw_rectangle(pmin=(start_rect, 0), pmax=(cur_pos, CHAR_HEIGHT), color=WHITE, fill=WHITE)
                            start_rect=0 #Reset to make new rectangle
                    cur_pos+=1
        with dpg.tooltip(f"{fasta_file}_{seq.id}_zoom_{cluster_changes_count}"):
            dpg.add_text(f"{seq.id}")
        

    def zoom_out_fasta_cluster(self, sender, app_data, user_data):
        """
        Shows a zoom out version of the aligned cluster fasta files represented with a line indicating the length of the
        sequence and rectangles indicating where alignments occur (or the pressence of a base in the fasta file)

        Input:
            sender = field dearpygui uses to bind item using the callback to the function zoom_out_fasta()
            app_data = field dearpygui uses
            user_data = file name at index 0; path to file at index 1
        
        Return:
            None
        """

        alias = user_data[0]
        zoom_alias = alias+"_zoom"
        if not dpg.does_alias_exist(zoom_alias):
            dpg.configure_item(user_data[0], show=False)
            fasta_file = user_data[1]
            
            seqs = list(SeqIO.parse(fasta_file, "fasta"))
            # maxIDLength = max(len(seq.id) for seq in seqs)
            num_seqs = len(seqs)
            if num_seqs > 15:
                num_seqs = 15
            with dpg.window(label=zoom_alias, tag=zoom_alias, show=True, horizontal_scrollbar=True, width=1500, height=(46 * num_seqs), pos=(250, 250)):
                # Display reference sequence 
                ref_seq = seqs[0]
                self.zoom_out_simplified_seq(alias, ref_seq, cur_row=0)
                cur_row = 1
                for seq in seqs[1:]: # Skip reference sequence
                    with dpg.item_handler_registry(tag=f"{fasta_file}_{seq.id}_popup_{cluster_changes_count}"):
                        dpg.add_item_clicked_handler(callback=self.change_order_popup, user_data=(seq, alias, fasta_file)) #Open a popup window when clicked on to change the cluster number of the sequence
                    self.zoom_out_simplified_seq(fasta_file, seq, cur_row) # Make the sequence line
                    dpg.bind_item_handler_registry(f"{fasta_file}_{seq.id}_zoom_{cluster_changes_count}", f"{fasta_file}_{seq.id}_popup_{cluster_changes_count}")
                    cur_row+=1
        else:
            dpg.configure_item(zoom_alias, show=True)
            dpg.configure_item(alias, show=False)
            dpg.focus_item(zoom_alias)


    def display_cluster_fasta(self, sender, app_data, user_data):
        """
        Displays the clsuter files in the clustered folder

        Input: 
            sender = field dearpygui uses to bind item using the callback to the function display_file()
            app_data = field dearpygui uses 
            user_data = cluster file name at index 0; cluster fasta file path at index 1

        Return:
            None
        """

        # Maybe in the future change it to where it doesn't create brand new item handler registries and tags and instead reuses old aliases
        global cluster_changes_count
        cluster_changes_count += 1

        alias = user_data[0]
        if not dpg.does_alias_exist(alias): # Make window if not already
            # Read fasta file
            fasta_file = user_data[1]
            seqs = list(SeqIO.parse(fasta_file, "fasta"))

            # Make maximum window height in case fasta file has many sequences
            num_seqs = len(seqs)
            if num_seqs > 10:
                num_seqs = 10
            
            ref_seq = seqs[0] # RefSeq should be the first element of the seqs list

            # Used for aligning sequences up
            maxIDLength = max([len(seq.id) for seq in seqs])
        
            # Display the clusters with the ability for the user to change which cluster group the sequence belongs to 
            with dpg.window(label=alias, tag=alias, show=True, horizontal_scrollbar=True, width=1500, height=(46 * num_seqs), pos=(250, 250)):

                # Submit changes in cluster group changes and update the order or view the sequences in a simplified manner
                with dpg.menu_bar():
                    dpg.add_menu_item(label="Save changes", callback=self.update_order, user_data=(alias, fasta_file))
                    dpg.add_menu_item(label="Zoom Out", callback=self.zoom_out_fasta_cluster, user_data=(alias, fasta_file))
                
                # Display the reference sequence
                with dpg.item_handler_registry(tag=f"{ref_seq}_{fasta_file}_location_{cluster_changes_count}"):
                    dpg.add_item_clicked_handler(callback=self.copy_seq, user_data = str(ref_seq.seq))
                with dpg.group(horizontal=True, horizontal_spacing=0):
                    ref_seq_id = f"{ref_seq.id}:"
                    dpg.add_text(f"{ref_seq_id:<{maxIDLength+5}}")
                    with dpg.group(horizontal=True, horizontal_spacing=0, tag=f"{ref_seq.id}_{fasta_file}_seq_{cluster_changes_count}"):
                        self.display_alignment(ref_seq)
                    with dpg.tooltip(f"{ref_seq.id}_{fasta_file}_seq_{cluster_changes_count}"):
                        dpg.add_text(f"Click to copy ({ref_seq.id})")
                    dpg.bind_item_handler_registry(f"{ref_seq.id}_{fasta_file}_seq_{cluster_changes_count}", f"{ref_seq.id}_{fasta_file}_location_{cluster_changes_count}")
                
                # Display the rest of the sequences
                with dpg.child_window(width=(7*(maxIDLength+len(ref_seq.seq))+61)):
                    for seq in seqs[1:]: # Place rest of the sequences after ref_seq
                        cluster_pattern_group = cluster_pattern.findall(seq.id)[0] # Skips the list part of re.findall to get the only tuple inside of the list
                        seq_cluster_num = cluster_pattern_group[0]
                        seq_name = cluster_pattern_group[1]

                        with dpg.group(horizontal=True, horizontal_spacing=0): #TODO: Figure out how to push input text and 'C' character to be right aligned instead of leaving large space
                            dpg.add_text("C")
                            dpg.add_input_text(default_value=seq_cluster_num, tag=seq_name, width=20, decimal=True) # Allow for user to change cluster group number
                            seq_id = f"{seq_name}:"
                            dpg.add_text(f"{seq_id:<{maxIDLength}}")
                            with dpg.item_handler_registry(tag=f"{seq_name}_{fasta_file}_location_{cluster_changes_count}"):
                                dpg.add_item_clicked_handler(callback=self.copy_seq, user_data=str(seq.seq))
                            with dpg.group(horizontal=True, horizontal_spacing=0, tag=f"{seq_name}_{fasta_file}_seq_{cluster_changes_count}"):
                                self.display_alignment(seq)
                            with dpg.tooltip(f"{seq_name}_{fasta_file}_seq_{cluster_changes_count}"):
                                dpg.add_text(f"Click to copy ({seq.id})")
                            dpg.bind_item_handler_registry(f"{seq_name}_{fasta_file}_seq_{cluster_changes_count}", f"{seq_name}_{fasta_file}_location_{cluster_changes_count}")

        else: # Otherwise display the window or focus it
            dpg.configure_item(alias, show=True)
            dpg.focus_item(alias)


    def display_file_txt(self, sender, app_data, user_data):
        """
        Displays the txt file

        Input:
            sender = field dearpygui uses to bind item (button) using the callback to the function display_file_txt
            app_data = field dearpygui uses
            user_data = tuple containing name of the cluster file at index 0 and path to the file as a string to index 1

        Return:
            None

        """
        
        alias = user_data[0]
        if not dpg.does_alias_exist(alias):
            # Get the number of lines in the txt file
            fasta_file = open(user_data[1])
            num_lines = len(fasta_file.readlines())
            fasta_file.close()
            if num_lines > 15:
                num_lines = 15

            #Display the content of the txt file
            fasta_file = open(user_data[1])
            with dpg.window(label=alias, tag=alias, show=True, horizontal_scrollbar=True, width=1080, height=num_lines*30, pos=(250, 250)):
                dpg.add_text(fasta_file.read())
            fasta_file.close()
        else:
            dpg.configure_item(alias, show=True)
            dpg.focus_item(alias)


    def display_file_fasta_aligned_specimens(self, sender, app_data, user_data):
        """
        Displasy the aligned fasta files for host and reference pair

        Input: 
            sender = field dearpygui uses to bind item (button) using the callback to the function display_file_fasta_aligned_sequences()
            app_data = field dearpygui uses
            user_data = tuple containing name of fasta file at index 0 and path to file as a string at index 1
        
        Return:
            None
        """

        alias = user_data[0]
        if not dpg.does_alias_exist(alias): # Make window if not already
            # Parse the fasta file
            fasta_file = f"{user_data[1]}"
            seqs = list(SeqIO.parse(fasta_file, "fasta"))
            maxIDLength = max(len(seq.id) for seq in seqs)
            
            # Create the window
            num_seqs = len(seqs)
            if num_seqs > 10:
                num_seqs = 10
            counter = 0
            with dpg.window(label=alias, tag=alias, show=True, horizontal_scrollbar=True, width=1500, height=30 + 25*num_seqs, pos=(250, 250)):
                # dpg.bind_font(courier_font)
                # Use a table to allow for readjustable vertical bar
                is_ref = False
                for seq in seqs:
                    with dpg.group(horizontal=True):
                        dpg.add_button(label=f"{seq.id:>{maxIDLength}}:", callback=self.copy_seq, user_data=str(seq.seq))
                        with dpg.tooltip(dpg.last_item()):
                            dpg.add_text("Click to copy to clipboard")

                        with dpg.item_handler_registry(tag=f"{seq.id}_{fasta_file}_location_{counter}"):
                            dpg.add_item_clicked_handler(button=dpg.mvMouseButton_Left, callback=self.copy_seq, user_data=str(seq.seq))
                        with dpg.group(horizontal=True, horizontal_spacing=0, tag=f"{seq.id}_{fasta_file}_seq_{counter}"): # Make the bases shown be horizontal with each other
                            self.display_alignment(seq)
                        with dpg.tooltip(parent=f"{seq.id}_{fasta_file}_seq_{counter}"):
                            dpg.add_text(f"Click to copy ({seq.id})")
                        dpg.bind_item_handler_registry(f"{seq.id}_{fasta_file}_seq_{counter}", f"{seq.id}_{fasta_file}_location_{counter}")
                    counter += 1
                    # Create separators between pairs of reference and host
                    if is_ref:
                        dpg.add_separator()
                        is_ref = False
                    else:
                        is_ref = True
        else: # Otherwise display the window
            dpg.configure_item(alias, show=True)
            dpg.focus_item(alias)


    def display_specimens_dir(self):
        """
            Helper function to display the specimens directory
        """

        specimens = self.file_sys_dictionary[RESULTS_DIR[:-1]][SPECIMEN_DIR[:-1]] # Omit the '/'
        for specimen in specimens:
            specimen_dir = specimens[specimen]
            with dpg.tree_node(label=specimen):
                specimen_path = f"{SPECIMEN_RESULTS_DIR}{specimen}/"
                with dpg.tree_node(label=DIAGRAMS_DIR[:-1]):
                    diagrams_dir = specimen_dir[DIAGRAMS_DIR[:-1]]
                    for v_family in diagrams_dir:
                        with dpg.tree_node(label=v_family):
                            for contig in diagrams_dir[v_family]: 
                                dpg.add_button(label=contig, callback=self.display_file_pdf, user_data=(contig, f"{specimen_path}{DIAGRAMS_DIR}{v_family}/{contig}"))
                with dpg.tree_node(label=SEQUENCES_DIR[:-1]):
                    sequences = specimen_dir[SEQUENCES_DIR[:-1]]
                    for sequence in sequences:
                        if "_aligned" in sequence:
                            dpg.add_button(label=sequence, callback=self.display_file_fasta_aligned_specimens, user_data=(sequence, f"{specimen_path}{SEQUENCES_DIR}{sequence}"))
                        else:
                            dpg.add_button(label=sequence, callback=self.display_file_fasta_not_aligned, user_data=(sequence, f"{specimen_path}{SEQUENCES_DIR}{sequence}"))


    def display_viruses_dir(self):
        """
            Helper function to display the viruses directory
        """

        viruses_dir = self.file_sys_dictionary[RESULTS_DIR[:-1]][VIRUS_DIR[:-1]] # Omit the '/'

        with dpg.tree_node(label=DIAGRAMS_DIR[:-1]):
            diagrams_dir = viruses_dir[DIAGRAMS_DIR[:-1]]
            for v_family in diagrams_dir.keys():
                with dpg.tree_node(label=v_family):
                    for virus in diagrams_dir[v_family].keys():
                        dpg.add_button(label=virus, callback=self.display_file_pdf, user_data=(virus, f"{v_diagrams_path}{v_family}/{virus}"))

        with dpg.tree_node(label=SEQUENCES_DIR[:-1]): 
            sequences_dir = viruses_dir[SEQUENCES_DIR[:-1]]
            for v_family in sequences_dir.keys():
                v_family_path = v_sequences_path + v_family
                with dpg.tree_node(label=v_family):
                    for virus in sequences_dir[v_family].keys():
                        if virus=="clustered": 
                            with dpg.tree_node(label=virus): #Tree node labeled clustered
                                for cluster_group in sequences_dir[v_family][virus].keys():
                                    cluster_group_file_path = f"{v_family_path}/{virus}/{cluster_group}"
                                    if sequences_dir[v_family][virus][cluster_group] != None:
                                        with dpg.tree_node(label=cluster_group):
                                            for cluster in sequences_dir[v_family][virus][cluster_group]:
                                                if fasta_pattern.fullmatch(cluster): # If the file is a fasta file
                                                    dpg.add_button(label=cluster, callback=self.display_file_fasta_not_aligned, user_data=(cluster, f"{cluster_group_file_path}/{cluster}"))
                                                elif txt_pattern.fullmatch(cluster): # If the file is a txt file
                                                    dpg.add_button(label=cluster, callback=self.display_file_txt, user_data=(cluster, f"{cluster_group_file_path}/{cluster}"))
                                                elif pdf_pattern.fullmatch(cluster): # If the file is a pdf file
                                                    dpg.add_button(label=cluster, callback=self.display_file_pdf, user_data=(cluster, f"{cluster_group_file_path}/{cluster}"))
                                    else:
                                        if fasta_pattern.fullmatch(cluster_group_file_path): # Determine if file is a fasta file
                                            dpg.add_button(label=cluster_group, tag=f"{cluster_group}_button", callback=self.display_cluster_fasta, user_data=(cluster_group, cluster_group_file_path))
                                        else: # Otherwise display the txt file version
                                            dpg.add_button(label=cluster_group, callback=self.display_file_txt, user_data=(cluster_group, cluster_group_file_path))
                        else:
                            if "_aligned" in virus:
                                dpg.add_button(label=virus, callback=self.display_file_fasta_aligned_viruses, user_data=(virus, f"{v_family_path}/{virus}"))
                            else:
                                dpg.add_button(label=virus, callback=self.display_file_fasta_not_aligned, user_data=(virus, f"{v_family_path}/{virus}"))


    def show_window(self):
        """
            Displays the results of EVE-X
        """

        dpg.create_context()
        dpg.create_viewport(title="EXE-X Results", width=2000, height=800)
        dpg.set_global_font_scale(1.0)
        # with dpg.font_registry():
        #     font = dpg.add_font("courier-normal.ttf", 13)
        with dpg.theme() as global_theme:
            with dpg.theme_component(dpg.mvAll):
                dpg.add_theme_style(dpg.mvStyleVar_ChildBorderSize, 0, category=dpg.mvThemeCat_Core)
        
        with dpg.window(tag="EVE-X Results", height=800):
            # dpg.bind_font(font)
            with dpg.table(header_row=False, resizable=True, height=800):
                dpg.add_table_column(init_width_or_weight=0.3)
                dpg.add_table_column()
                with dpg.table_row():
                    with dpg.table_cell():
                        with dpg.child_window(horizontal_scrollbar=True):
                            with dpg.collapsing_header(label=SPECIMEN_DIR[:-1]):
                                self.display_specimens_dir()
                            with dpg.collapsing_header(label=VIRUS_DIR[:-1]):
                                self.display_viruses_dir()
        dpg.setup_dearpygui()
        # dpg.show_debug()
        # dpg.show_item_registry()
        dpg.set_exit_callback(callback=self.exit_app)
        dpg.show_viewport()
        dpg.set_primary_window("EVE-X Results", True)
        dpg.start_dearpygui()
        dpg.destroy_context()


def main():
    # Verify Paths are spelled correctly and that the directories exist
    p = Path(SPECIMEN_RESULTS_DIR)
    if not p.exists():
        sys.exit("Error: Specimens Results Directory does not exist")
    p = Path(VIRUS_RESULTS_DIR)
    if not p.exists():
        sys.exit("Error: Virus Results Directory does not exist")
    
    p = Path(RESULTS_DIR)
    Interface(p)

if __name__ == "__main__":
    main()
