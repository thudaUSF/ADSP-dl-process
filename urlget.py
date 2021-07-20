import sys
import boto3

snd = sys.argv[1]
ID = sys.argv[2]

client = boto3.client('s3')
url = client.generate_presigned_url("get_object", Params={"Bucket":"wanglab-dss-share","Key":'distribution/adsp/cram/' + str(snd) + '/' + str(ID) + '.cram', "RequestPayer":'requester'}) #changes depending on which cram you download
print("./htsfile -h '{0}'".format(url)) 
sys.exit(0)
