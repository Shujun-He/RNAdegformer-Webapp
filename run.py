# Form / Frame
# Use a frame component in a form card to display HTML content inline.
# ---
import sys
sys.path.append('draw_rna_pkg/')
from h2o_wave import Q, listen, ui
import matplotlib.pyplot as plt
from ipynb.draw import draw_struct
import os
import io
import base64
import matplotlib.image as mpimg
import zipfile
import time



os.environ["ARNIEFILE"] = f"arnie.conf"


from rna_analysis import *
from dna_analysis import *
import RNA_Inference

#from model_prediction import load_models,bpps_path,get_prediction_df_dict,bpps_check






# #for cmaps check out https://matplotlib.org/3.2.1/tutorials/colors/colormaps.html
target_columns = [ 'reactivity', 'deg_Mg_pH10',
       'deg_pH10', 'deg_Mg_50C', 'deg_50C']

logo_file = "files/logo.png"
image_1_path = "files/RNAdegformer.png"
image_2_path = "files/image_2_path.png"
image_3_path = "files/image_3_path.png"
image_4_path = "files/image_4_path.png"

all_pages = ['nav','home', 'upload', 'example', 'infoboard1', 'image0', 'image1', 'image2','data_info','stn_hist',
                 'image3', 'image4', 'image5', 'description', 'text','text2', 'error', 'progress', 'selection_card', 'plot',
                 "plot_seqpos","plot_count",'plot_features0','plot_features1','plot_features2','plot_features3','text_sub','predict_tab','predict_promoter_tab']


os.system('mkdir temp')

image_size = 6
number_of_plots_in_a_row = 3
plot_height = 5
plot_width = 3

data_display_max_nrows = 10



def get_image(file_name = "image1"):
    plt.figure(figsize=(image_size, image_size))
    buf = io.BytesIO()

    img = mpimg.imread(f'temp/{file_name}.png')
    plt.imshow(img)
    plt.axis('off')
    plt.savefig(buf, format='png')
    buf.seek(0)
    image = base64.b64encode(buf.read()).decode('utf-8')

    return image

def display_error(error_type):
    items = []
    if error_type == "upload":
        items = [ui.message_bar(type='error', text='**Selected file was not in the expected format.**'),
                    ui.button(name='upload_data', label='Upload again', primary=True),
                    ui.separator()]
    elif error_type == "misentered_info":
        items = [ui.message_bar(type='error', text='**Either sequence or structure info not entered.**'),
                ui.button(name='custom_back', label='Try again', primary=True),
                ui.separator()]
    elif error_type == "wrong_columns":
        items = [
                ui.message_bar(type='error', text='**Either sequence or structure column not selected properly**'),
                ui.button(name='restart', label='Try again', primary=True),
                ui.separator()
            ]
    elif error_type == "misentered_ids":
        items = [
                ui.message_bar(type='error', text='**IDs not entered properly**'),
                ui.button(name='restart', label='Try again', primary=True),
                ui.separator()
            ]
    elif error_type == "arnie_not_selected":
        items = [
            ui.message_bar(type='error', text='**Either feature or package not selected properly.**'),
            ui.button(name='arnie_back', label='Try again', primary=True),
            ui.separator()
        ]
    return items

def make_ui_table(df, n_rows, name = "head_of_table"):
    """Creates a ui.table object from a csv file"""

    n_rows = min(n_rows, df.shape[0])

    table = ui.table(
        name=name,
        columns=[ui.table_column(name=str(x), label=str(x), sortable=True) for x in df.columns.values],
        rows=[ui.table_row(name=str(i), cells=[str(df[col].values[i]) for col in df.columns.values])
              for i in range(n_rows)]
    )
    return table

async def delete_pages(q: Q,keep_nav=False):
    # There does not seem to be a way to know what pages are in q.page
    page_list = ['nav','home', 'upload', 'example', 'infoboard1','data_info','stn_hist',
                 'description', 'text','text2', 'error', 'progress', 'selection_card', 'plot','data_view',
                 "plot_seqpos","plot_count",'plot_features0','plot_features1','plot_features2','plot_features3','text_sub','predict_tab','predict_promoter_tab'] \
                + list(set(q.client.all_pages))
    if keep_nav:
        page_list.remove("nav")

    for page in page_list:
        try:
            del q.page[page]
        except:
            #print('some error')
            pass

async def display_nav(q: Q):



    q.page['nav'] = ui.tab_card(
        box='3 1 9 1',
        items=[
            ui.tab(name='#home', label='Home'),
            ui.tab(name='#rnaprediction', label='RNA Degradation Prediction'),
        ],
        link=False
    )

    q.page['header2'] = ui.header_card(
        box='1 1 2 1',
        title='RNAdegformer',
        subtitle='',
        icon='ExploreData',
        icon_color='$orange',
    )

    q.page['logo'] = ui.markup_card(
        box='12 1 1 1',
        title='',
        content="""<p style='text-align:center; vertical-align: middle; display: table-cell; width: 134px;'>"""
                """<a href='https://www.h2o.ai/products/h2o-wave/'> <img src='""" + q.app.logo_url + """' height='50px' width='50px'> </a> </p>""" )

async def progress_page(q: Q, message='\nJust a second...'):
    q.page['progress'] = ui.form_card(box='4 5 6 2',
                                      items=[ui.progress(label=message)])
    await q.page.save()
    del q.page['progress']



async def display_file_upload(q):
    print("location: display_file_upload")

    q.page['description'] = ui.form_card(box='1 2 3 10',
        items=[ui.text_m('\n You can upload a local dataset to run ExploRNA.'),
               ui.file_upload(name='user_files', label='Upload', multiple=False)])

    if q.client.train is not None:
        data_items = [ui.text_m(f'Loaded file "{q.client.file_name}" has '
                                f'**{q.client.train.shape[0]}** rows and **{q.client.train.shape[1]}** features.\n\n'),
                      make_ui_table(q.client.train, data_display_max_nrows)]
        q.page['data_view'] = ui.form_card(box='4 2 9 7', items=data_items)

















async def predict_rna_tab(q):
    ui_list_custom = [
        ui.text_m(f'In this section, you can directly type a sequence and get predictions for the following targets: {target_columns}.'),
        ui.text_m(f'Additional features will be generated automatically and we will visulize the RNA folding and \
        attention weights of the Nucleic Transformer for you as well'),
        ui.textbox(name='rna_sequence_textbox', label='Sequence', value=q.args.rna_sequence_textbox or "GGAAAAGCUCUAAUAACAGGAGACUAGGACUACGUAUUUCUAGGUAACUGGAAUAACCCAUACCAGCAGUUAGAGUUCGCUCUAACAAAAGAAACAACAACAACAAC"),
        ui.button(name='predict_rna', label='Predict', primary=True),
        ui.text_s(f'You can enter in following formats: '),
        ui.text_s(f'Sequence characters: A,U,C,G'),
    ]
    q.page['text'] = ui.form_card( box='1 2 3 4',items=ui_list_custom )

    if q.args.rna_sequence_textbox is not None:
        q.client.rna_sequence_textbox=q.args.rna_sequence_textbox

    if q.client.rna_predictions is not None and q.client.rna_sequence_textbox is not None:

        # q.client.rna_drawn_html_pos = position_based_plot_single(q.client.rna_predictions, target_columns, size=500)
        # q.page['plot_seqpos'] = ui.frame_card(
        #     box='4 2 9 4',
        #     title='Average prediction value plots for each sequence position.',
        #     content=q.client.rna_predictions_html_pos
        # )



        q.client.rna_predictions_html_pos = position_based_plot_single(q.client.rna_predictions, target_columns, size=500)
        q.page['plot_seqpos'] = ui.frame_card(
            box='4 7 9 6',
            title='Average prediction value plots for each sequence position.',
            content=q.client.rna_predictions_html_pos
        )

        draw_struct(q.client.rna_sequence_textbox, q.client.rna_input_features['structures'][0], c=None, ax=None,file_name='rna_plot')
        image = get_image(file_name='rna_plot')

        #image_path_list = [f"temp/{custom_id}.png"]
        #compress(image_path_list)
        # Deleting previously generated plots after compressing them.
        # for path in image_path_list:
        #     os.remove(path)

        #q.client.plots_path, = await q.site.upload([f'temp/plots.zip'])

        q.page['image1'] = ui.image_card(
            box=f'4 2 4 5',
            title=f'RNA Visualization',
            type='png',
            image=image)

        q.client.rna_aw_html = plot_aw(q.client.rna_aw, size=500)

        q.page['aw_image'] = ui.frame_card(
            box=f'8 2 5 5',
            title=f'Attetion weight of the Nucleic Transformer',
            content=q.client.rna_aw_html)

        q.client.all_pages.append("image1")
        q.client.all_pages.append("aw_image")

        download_data_text = '''=
        Inference complete! Click [here]({{predictions}}) to download the predictions!
        '''
        # q.page['download_rna_predictions'] = ui.markdown_card(
        #         box='8 2 5 5',
        #         title='',
        #         content=download_data_text,
        #         data=dict(predictions=q.client.rna_predictions_path)
        #     )

        download_data_text = '''=
Inference complete! Click [here]({{predictions}}) to download the predictions!
'''
        q.page['download_promoter_predictions'] = ui.markdown_card(
                box='1 6 3 1',
                title='',
                content=download_data_text,
                data=dict(predictions=q.client.rna_predictions_path)
            )


async def rna_model_predict(q):

    # Loading models
    start_time_main = time.time()

    if q.client.rna_models_loaded is None:
        q.page['progress'] = ui.form_card(box='4 6 9 1',
                                          items=[ui.progress(label="Loading models...")])
        await q.page.save()
        del q.page['progress']
        q.client.rna_inference=RNA_Inference.RNA_Inference()
        q.client.rna_inference.load_models('RNA_Inference/best_weights')
        q.client.rna_models_loaded=True

    # Getting predictions
    #if q.client.predictions is None:
    #minimum_required_features = ["sequence"]
    #await progress_page(q,message="Getting predictions...")

    q.page['progress'] = ui.form_card(box='4 6 9 1',
                                      items=[ui.progress(label="Getting predictions...")])
    await q.page.save()
    del q.page['progress']

    q.client.rna_predictions, q.client.rna_aw, q.client.rna_input_features  = \
    q.client.rna_inference.predict(q.args.rna_sequence_textbox)
    q.client.rna_predictions_df=pd.DataFrame(columns=['position']+target_columns)
    q.client.rna_predictions_df['position']=np.arange(len(q.client.rna_predictions))
    q.client.rna_predictions_df[target_columns]=q.client.rna_predictions
    q.client.rna_predictions_df.to_csv('temp/rna_predictions.csv', index=False)
    q.client.rna_predictions_path, = await q.site.upload(['temp/rna_predictions.csv'])

    #plt.imshow(q.client.rna_aw)
    #plt.savefig('temp/aw.png',bbox_inches = 'tight',pad_inches = 0)

    elapsed_time = time.time() - start_time_main
    print(f"minutes passed for getting predictions: {round(elapsed_time/60, 2)}")


async def home(q):
    home_text = f'''

The RNAdegformer is a deep learning model developed to study and understand RNA degradation. The model archiecture is simple but effective, and
we used it to to place 7th during the [OpenVaccine challenge](https://www.kaggle.com/c/stanford-covid-vaccine) and outperform the top solution with a semi-supervised version later.

Throughout the app, you will be able to use trained RNAdegformer models to quantify RNA degradation and visualize their secondary structure. Also,
Nucleic Transformer is very interpretible so you will be able to visualize the attention of the neural networks and understand what the neural network is looking at/for while making predictions.
This also allows informed decision making when using these neural networks.

![Plot]({q.app.home_image_1_url})

 '''
    q.page['home'] = ui.form_card(box='1 2 12 10',
        items=[
               ui.text_m(home_text)
               ])

async def display_data_parameters_page(q,random_sample_disabled=True):

    # Add navigation buttons
    print(q.args["#"])

    if q.client.activetab == "home":
            await delete_pages(q,keep_nav=True)
            await home(q)
    elif q.client.activetab == "rnaprediction":
            await delete_pages(q,keep_nav=True)
            await predict_rna_tab(q)






async def main(q: Q):
    # Upload the logo & images
    if q.client.all_pages is None:
        q.client.all_pages = []

    if not q.app.logo_url:
        q.app.logo_url, = await q.site.upload([logo_file])

    if not q.app.home_image_1:
        q.app.home_image_1_url, = await q.site.upload([image_1_path])

    if not q.app.home_image_2:
        q.app.home_image_2_url, = await q.site.upload([image_2_path])

    if not q.app.home_image_3:
        q.app.home_image_3_url, = await q.site.upload([image_3_path])

    if not q.app.home_image_4:
        q.app.home_image_4_url, = await q.site.upload([image_4_path])


    q.client.activetab = 'home'
    await display_data_parameters_page(q)
    await display_nav(q)


    if q.args.predict_rna:
        #q.client.virus_topk = q.args.virus_topk
        await rna_model_predict(q)
        print('line 815')
        await predict_rna_tab(q)


    elif q.args["#"]:
        q.client.activetab = q.args["#"]
        await display_data_parameters_page(q)

    else:
        print("location: else")
        await display_data_parameters_page(q)

    await q.page.save()


if __name__ == '__main__':
    listen('/rnadegformer', main)
