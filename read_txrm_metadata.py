import olefile
import struct
import sys
import numpy as np
from libtiff import TIFF
import os
from datetime import datetime

def read_txrm_metadata(in_filename,basepath,basefilename):
    ole=olefile.OleFileIO(in_filename);
    
    if (ole.exists('ImageInfo/NoOfImages') and 
                ole.exists('ImageInfo/ImageWidth') and 
    ole.exists('ImageInfo/ImageHeight')):                  
        stream = ole.openstream('ImageInfo/NoOfImages')
        data = stream.read()
        nimages = struct.unpack('<I', data)
        NImages=nimages[0]

        stream = ole.openstream('ImageInfo/ImageHeight')
        data = stream.read()
        ximage = struct.unpack('<I', data)    
        Height= np.int(ximage[0]);
        
        stream = ole.openstream('ImageInfo/ImageWidth')
        data = stream.read()
        yimage = struct.unpack('<I', data)
        Width= np.int(yimage[0]);

        stream = ole.openstream( 'ImageInfo/pixelsize')
        data = stream.read()
        pixelsize = struct.unpack('<f', data)
        PixelSize= np.float(pixelsize[0]);

#        stream = ole.openstream( 'ImageInfo/ExpTimes')
#        data = stream.read()
#        ExpTimes = struct.unpack('<f', data)
#        ExpTimes= np.float(ExpTimes[0]);
#        
#        stream = ole.openstream( 'ImageInfo/Voltage')
#        data = stream.read()
#        Voltage = struct.unpack('<f', data)
#        Voltage= np.float(Voltage[0]);
#
#        stream = ole.openstream( 'ImageInfo/RequestedPower')
#        data = stream.read()
#        RequestedPower = struct.unpack('<f', data)
#        RequestedPower= np.float(RequestedPower[0]);
#
#        stream = ole.openstream( 'ImageInfo/SourceFilterName')
#        data = stream.read()
#        SourceFilterName = struct.unpack('<s', data)
#        SourceFilterName= np.str(SourceFilterName[0]);


    DataType='float'
    if ole.exists('ImageInfo/DataType'):                  
        stream = ole.openstream('ImageInfo/DataType')
        data = stream.read()
        struct_fmt = '<1I'
        datatype = struct.unpack(struct_fmt, data)
        datatype = int(datatype[0])
        if datatype == 5:
            DataType = 'uint16'
        else:
            DataType = 'float'

    if ole.exists('ImageInfo/Date'):   
        stream = ole.openstream('ImageInfo/Date')       
        data = bytearray(stream.read());
#        dates = struct.unpack('<'+'17s23x'*NImages, data) 
        dates = struct.unpack('<'+'19s21x'*NImages, data)
        for iimage in range(NImages):
        
            #Read the images - They are stored in the txrm as ImageData1, ImageData2... 
            #Each folder contains 100 images 1-100, 101-200... 
            img_string = "ImageData%i/Image%i" % (np.ceil((iimage+1)/100.0), (iimage+1))
            print img_string
            stream = ole.openstream(img_string)
            data = stream.read()
        
            if DataType == 'uint16':
                struct_fmt = "<{0:10}H".format(Height*Width)
                imgdata = struct.unpack(struct_fmt, data)
                imgdata=np.uint16(imgdata)
            elif DataType == 'float':                   
                struct_fmt = "<{0:10}f".format(Height*Width)
                imgdata = struct.unpack(struct_fmt, data)
            else:                            
                print "Wrong data type"
                
            singleimage = np.flipud(np.reshape(imgdata, (Height, Width), order='A'))
            singleimage = np.reshape(singleimage, (1, Height, Width), order='A')
            #singleimage = np.reshape(imgdata, (1, Height, Width), order='A')
            filename=dates[iimage];
            filename=str.replace(filename,'/','');
            filename=str.replace(filename,' ','-');
            filename=str.replace(filename,':','');
            try:
                tempdate=datetime.strptime(dates[iimage],'%m/%d/%Y %H:%M:%S') 
                tempfilename=os.path.join(basepath,basefilename + '_' +  \
                    datetime.strftime(tempdate,'%y%m%d%H%M%S'))
            except:
                tempfilename=os.path.join(basepath,basefilename + '_')
            imfilename=os.path.join(tempfilename+'.tif')   
            exists = os.path.isfile(imfilename)
            print "test {0}".format(exists) + ", {0}".format(imfilename)
            iversion=1
            while exists:
                try:
                    tempfilename=os.path.join(basepath,basefilename + '_' +  \
                        datetime.strftime(tempdate,'%y%m%d%H%M%S') + '_{0:0>4}'.format(iversion))
                except:
                    tempfilename=os.path.join(basepath,basefilename + '_{0:0>4}'.format(iversion))
                imfilename=os.path.join(tempfilename+'.tif')   
                exists = os.path.isfile(imfilename)
                print "test {0}".format(exists) + ", {0}".format(imfilename)
                iversion+=1
                
            metadatafile=os.path.join(tempfilename+'.txt')
            
            print "writing metadata: " + metadatafile
            file = open(metadatafile,'w') 
            tempstr= "ImageInfo>ImageHeight = %i" % Height 
            print tempstr
            file.write(tempstr+'\n')   
            tempstr= "ImageInfo>ImageWidth = %i" % Width
            print tempstr
            file.write(tempstr+'\n')      
            tempstr= "ImageInfo>pixelsize = %f" % PixelSize
            print tempstr
            file.write(tempstr+'\n')   
            tempstr= "ImageInfo>DataType = %s " % DataType 
            print tempstr
            file.write(tempstr+'\n')   
            tempstr= "ImageInfo>Date = %s " % dates[iimage]
            file.write(tempstr+'\n')   
            file.close() 
            
                        
            print "writing image: " + imfilename
            tiff = TIFF.open(imfilename, mode='w')
            tiff.write_image(singleimage)
            tiff.close()
        

if __name__ == '__main__':
    in_filename = str(sys.argv[1])
    basepath = str(sys.argv[2])
    basefilename = str(sys.argv[3])
    sys.stdout.write(str(read_txrm_metadata(in_filename,basepath,basefilename)))


        