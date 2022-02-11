def apply(st):
    # Hide the menu button
    st.markdown(""" <style>
                #MainMenu {visibility: hidden;}
                footer {visibility: hidden;}
                </style> """, unsafe_allow_html=True)
    # Do some css styling tricks here (e.g. remove the padding)
    # https://medium.com/ssense-tech/streamlit-tips-tricks-and-hacks-for-data-scientists-d928414e0c16
    padding = 1
    st.markdown(f""" <style>
                .css-12oz5g7{{
                padding-top: {padding}rem;
                margin-top: -3.5rem;
                max-width: 50rem;
                padding-right: {padding}rem;
                padding-left: {padding}rem;
                padding-bottom: {padding}rem;
                }} </style> """, unsafe_allow_html=True)
    st.markdown(f""" <style>
                .css-128j0gw{{
                margin-top: -3.0rem;
                }} </style> """, unsafe_allow_html=True)
    # Reduce the space in horizontal component
    st.markdown(f""" <style>
                hr {{
                margin: 10px 0px;
                }} </style> """, unsafe_allow_html=True)
    # Custom button appearance
    st.markdown('''
                <style>
                div.stButton > button:first-child {
                display: inline-flex;
                align-items: center;
                justify-content: center;
                min-width: 100%;
                background-color: rgb(255, 255, 255);
                color: rgb(38, 39, 48);
                padding: .25rem .75rem;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(230, 234, 241);
                border-image: initial; }
                div.stButton > button:hover {
                border-color: rgb(246, 51, 102);
                color: rgb(246, 51, 102);}
                div.stButton > button:active{
                box-shadow: none;
                background-color: rgb(246, 51, 102);
                color: white; }
                </style>''', unsafe_allow_html=True)
    # Custom button appearance
    st.markdown('''
                <style>
                div.stDownloadButton > button:first-child {
                display: inline-flex;
                align-items: center;
                justify-content: center;
                min-width: 100%;
                background-color: rgb(255, 255, 255);
                color: rgb(38, 39, 48);
                padding: .25rem .75rem;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(230, 234, 241);
                border-image: initial;}
                div.stDownloadButton > button:hover {
                border-color: rgb(246, 51, 102);
                color: rgb(246, 51, 102);}
                div.stDownloadButton > button:active {
                box-shadow: none;
                background-color: rgb(246, 51, 102);
                color: white;}
                </style>''', unsafe_allow_html=True)
