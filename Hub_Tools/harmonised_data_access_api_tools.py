import requests, re, json, urllib3, sys
import shutil
import time, os
import urllib.parse

def init(dataset_id, api_key, download_dir_path, verbose=True):
    '''
     Initisalise Harmonised Data Access API path dictionary
    '''
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    HAPI_dict = {}
    # HDA-API endpoint
    HAPI_dict["apis_endpoint"]="https://apis.wekeo.eu"
    # Data broker address
    HAPI_dict["broker_address"] = HAPI_dict["apis_endpoint"]\
                                  + "/databroker/0.1.0"
    # Terms and conditions
    HAPI_dict["acceptTandC_address"]\
            = HAPI_dict["apis_endpoint"]\
            + "/dcsi-tac/0.1.0/termsaccepted/Copernicus_General_License"
    # Access-token address
    HAPI_dict["accessToken_address"] = HAPI_dict["apis_endpoint"]\
                                       + '/token'
    # Dataset id
    HAPI_dict["dataset_id"] = dataset_id
    HAPI_dict["encoded_dataset_id"] = urllib.parse.quote(dataset_id)
    # API key
    HAPI_dict["api_key"] = api_key

    # set HTTP success code
    HAPI_dict["CONST_HTTP_SUCCESS_CODE"] = 200

    # download directory
    HAPI_dict["download_dir_path"] = download_dir_path
    if not os.path.exists(download_dir_path):
        os.makedirs(download_dir_path)

    # set verbosity
    HAPI_dict["verbose"] = verbose

    return HAPI_dict

def get_access_token(HAPI_dict):
    '''
     Get an access token from an API key
    '''
    headers = {'Authorization': 'Basic ' + HAPI_dict["api_key"]}
    data = [('grant_type', 'client_credentials'), ]
    print("Getting an access token. This token is valid for one hour only.")
    response = requests.post(HAPI_dict["accessToken_address"],\
                headers=headers, data=data, verify=False)

    # If the HTTP response code is 200 (i.e. success), then retrive the token from the response
    if (response.status_code == HAPI_dict["CONST_HTTP_SUCCESS_CODE"]):
        access_token = json.loads(response.text)['access_token']
        print("Success: Access token is " + access_token)
    else:
        print("Error: Unexpected response {}".format(response))

    HAPI_dict["access_token"] = access_token

    HAPI_dict["headers"] = \
       {'Authorization': 'Bearer ' + HAPI_dict["access_token"],}

    return HAPI_dict

def query_metadata(HAPI_dict):
    '''
     Check data product metadata
    '''
    response = requests.get(HAPI_dict["broker_address"] + '/querymetadata/'\
               + HAPI_dict["encoded_dataset_id"],\
               headers=HAPI_dict["headers"])
    print('Getting query metadata, URL Is '\
          + HAPI_dict["broker_address"] + '/querymetadata/'\
          + HAPI_dict["encoded_dataset_id"] + "?access_token="\
          + HAPI_dict["access_token"])
    print("************** Query Metadata for " + HAPI_dict["dataset_id"]\
         +" **************")

    if (response.status_code == HAPI_dict["CONST_HTTP_SUCCESS_CODE"]):
        parsedResponse = json.loads(response.text)
        print(json.dumps(parsedResponse, indent=4, sort_keys=True))
        print("*****************************************"\
                  + "*********************************")
    else:
        print("Error: Unexpected response {}".format(response))
        parsedResponse = "Error"

    HAPI_dict["parsedResponse"] = parsedResponse
    return HAPI_dict

def accept_TandC(HAPI_dict):
    '''
     Accept terms a conditions (could replace with tick-box wiki, but not automated)
    '''
    response = requests.get(HAPI_dict["acceptTandC_address"],\
               headers=HAPI_dict["headers"])
    isTandCAccepted = json.loads(response.text)['accepted']
    if isTandCAccepted is False:
        print("Accepting Terms and Conditions of Copernicus_General_License")
        response = requests.put(acceptTandC_address, headers=headers)
    else:
        print("Copernicus_General_License Terms and Conditions already accepted")

    HAPI_dict["isTandCaccepted"] = True
    return HAPI_dict

def launch_query(HAPI_dict, query):
    '''
     Launch query and get job ID
    '''
    response = requests.post(HAPI_dict["broker_address"] + '/datarequest',\
               headers=HAPI_dict["headers"],\
               json=query, verify=False)
    if (response.status_code == HAPI_dict["CONST_HTTP_SUCCESS_CODE"]):
        job_id = json.loads(response.text)['jobId']
        print ("Query successfully submitted. Job ID is " + job_id)
    else:
        print("Error: Unexpected response {}".format(response))
        job_id = None

    HAPI_dict["job_id"] = job_id
    return HAPI_dict

def check_job_status(HAPI_dict):
    '''
     Check status of HAPI job
    '''
    isComplete = False
    while not isComplete:
        response = requests.get(HAPI_dict["broker_address"]\
                   + '/datarequest/status/' + HAPI_dict["job_id"],\
                   headers=HAPI_dict["headers"])
        results = json.loads(response.text)['resultNumber']
        isComplete = json.loads(response.text)['complete']
        if isComplete:
            print("The Job " + HAPI_dict["job_id"]\
              + " has completed")
        else:
           print("The Job " + HAPI_dict["job_id"]\
              + " has not completed")
           # sleep for 2 seconds before checking the job status again
           time.sleep(2)

    numberOfResults = str(results)
    HAPI_dict['nresults'] = results
    print ("Total number of products/results :" + numberOfResults)
    
    return HAPI_dict

def get_results_list(HAPI_dict):
    '''
     Return query results
    '''
    response = requests.get(HAPI_dict["broker_address"]\
               + '/datarequest/jobs/' + HAPI_dict["job_id"] + '/result',\
               headers=HAPI_dict["headers"])
    results = json.loads(response.text)

    if HAPI_dict["verbose"] == True:
        print("************** Results *******************************")
        print(json.dumps(results, indent=4, sort_keys=True))
        print("*********************************************")

    HAPI_dict["results"] = results
    return HAPI_dict

def get_download_links(HAPI_dict):
    '''
     Get download links associated with results
    '''
    download_urls = []
    for result in HAPI_dict["results"]['content']:
        externalUri = result['externalUri']
        product_size = result['fileSize']/(1024*1024)
        product_name = result['fileName']
        download_url = HAPI_dict["broker_address"]\
                       + '/datarequest/result/' + HAPI_dict["job_id"]\
                       + '?externalUri='\
                       + urllib.parse.quote(externalUri)\
                       +"&access_token="+HAPI_dict["access_token"]
        if HAPI_dict["verbose"] == True:
            print("Download link for " + product_name\
                  + "(" + "{:.2f}".format(product_size) + " MB) :")
            print(download_url)
        download_urls.append(download_url)

    HAPI_dict["download_urls"] = download_urls
    return HAPI_dict

def downloadFile(url, headers, directory, total_length = 0) :
    """
    downloadFile(url, headers, directory)
    Download the file using streaming
    Parameters
    ----------
    url:
        The URL for the request
    headers:
        Headers (e.g. auth bearer token)
    directory:
        Directory for the download
    size:
        The size of the file
    """
    r = requests.get(url, headers=headers, stream=True)
    if r.status_code == 200:
        filename = os.path.join(directory,\
                   get_filename_from_cd(r.headers.get('content-disposition')))
        print("Downloading " + filename)
        with open(filename, 'wb') as f:
            start = time.clock()
            print("File size is: %8.2f MB" % (total_length/(1024*1024)))
            dl = 0
            for chunk in r.iter_content(64738):
                dl += len(chunk)
                f.write(chunk)
                if total_length is not None:
                    done = int(50 * dl / total_length)
                    str1 = '=' * done
                    str2 =' ' * (50-done)
                    str3 = (dl//(time.clock() - start))/(1024*1024)
                    print("\r[%s%s] %8.2f Mbps" % (str1, str2, str3), end='', flush=True)
                else:
                    if( dl % (1024) == 0 ):
                        str1 = dl / (1024 * 1024)
                        str2 = (dl//(time.clock() - start))/1024
                        print("[%8.2f] MB downloaded, %8.2f kbps" % (str1, str2))
            str1 = dl / (1024 * 1024)
            str2 = (dl//(time.clock() - start))/1024
            print("[%8.2f] MB downloaded, %8.2f kbps" % (str1, str2))
            return (time.clock() - start)

def get_filename_from_cd(cd):
    """
    get_filename_from_cd(cd)
    Get the filename from content disposition
    Parameters
    ----------
    cd : content disposition
        From https://www.w3.org/Protocols/rfc2616/rfc2616-sec19.html
        The Content-Disposition response-header field has been proposed 
        as a means for the origin server to suggest a default filename
        if the user requests that the content is saved to a file.
        This usage is derived from the definition of Content-Disposition in RFC 1806 [35].
    """
    if not cd:
        return None
    fname = re.findall('filename=(.+)', cd)
    if len(fname) == 0:
        return None
    return fname[0].replace("'","").replace('"',"")

def download_data(HAPI_dict, skip_existing=False): 
    '''
     Download the data
    '''
    counter = 0
    filenames = []
    for result in HAPI_dict["results"]['content']:
        externalUri = result['externalUri']
        product_size = result['fileSize']
        download_url = HAPI_dict["broker_address"]\
                       + '/datarequest/result/'\
                       + HAPI_dict["job_id"] + '?externalUri='\
                       + urllib.parse.quote(externalUri)

        r = requests.get(download_url, headers=HAPI_dict["headers"],\
                         stream=True)
        filename = os.path.join(HAPI_dict["download_dir_path"],\
                   get_filename_from_cd(r.headers.get('content-disposition')))
        filenames.append(filename)
        
        if skip_existing and os.path.exists(filename):
            print("Skipping " + os.path.basename(filename) + " as it exists already")
        else:
            time_elapsed = downloadFile(download_url, HAPI_dict["headers"],\
                          HAPI_dict["download_dir_path"], product_size)
            print("Download complete (took " + str(time_elapsed) + " seconds)")
            print("")
            
    HAPI_dict['filenames'] = filenames
    return HAPI_dict

def get_filenames(HAPI_dict):
    '''
     Get the filenames of the target file
    '''
    filenames = []
    for result in HAPI_dict["results"]['content']:
        externalUri = result['externalUri']
        product_size = result['fileSize']
        download_url = HAPI_dict["broker_address"]\
                       + '/datarequest/result/'\
                       + HAPI_dict["job_id"] + '?externalUri='\
                       + urllib.parse.quote(externalUri)

        r = requests.get(download_url, headers=HAPI_dict["headers"],\
                         stream=True)
        filename = os.path.join(HAPI_dict["download_dir_path"],\
                   get_filename_from_cd(r.headers.get('content-disposition')))
        filenames.append(filename)
    HAPI_dict['filenames'] = filenames
    return HAPI_dict
    
    