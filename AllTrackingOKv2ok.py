import sys
import os
import ij
from ij import IJ, ImagePlus
from ij.io import DirectoryChooser
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from fiji.plugin.trackmate import Model, Settings, TrackMate, SelectionModel, Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.gui.displaysettings.DisplaySettings import TrackMateObject
from fiji.plugin.trackmate.features.track import TrackIndexAnalyzer
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.action import ExportTracksToXML
from fiji.plugin.trackmate.io import TmXmlWriter
from java.io import File
import java.io.FileWriter as FileWriter
import java.io.IOException as IOException

# Avoid errors with UTF8 chars in TrackMate
reload(sys)
sys.setdefaultencoding('utf-8')

# Ask the user to select the directory containing the .nd files
dc = DirectoryChooser("Select the directory containing .nd files")
input_dir = dc.getDirectory()

# Function to process a single .nd file
def process_nd_file(nd_file_path):
    print("Processing file: " + nd_file_path)
    
    # Get the base name of the input .nd file
    base_name = os.path.splitext(os.path.basename(nd_file_path))[0]

    # Parameters for each channel
    channel_params = [
        {'RADIUS': 3.5, 'THRESHOLD': 1.0, 'OUTPUT_XML': os.path.join(input_dir, '{}_w1Quad-TIRF488_t1.xml'.format(base_name)), 'OUTPUT_CSV': os.path.join(input_dir, '{}_channel_1_results.csv'.format(base_name))},
        {'RADIUS': 3.5, 'THRESHOLD': 1.0, 'OUTPUT_XML': os.path.join(input_dir, '{}_w2Quad-TIRF642_t1.xml'.format(base_name)), 'OUTPUT_CSV': os.path.join(input_dir, '{}_channel_2_results.csv'.format(base_name))}
    ]

    # 1. Open the .nd file using Bio-Formats
    options = ImporterOptions()
    options.setId(nd_file_path)
    options.setSplitChannels(True)
    imps = BF.openImagePlus(options)

    # Process each channel independently
    for i, imp in enumerate(imps):
        imp.show()

        # Enhance contrast for the image
        IJ.run(imp, "Enhance Contrast", "saturated=0.35")

        # Create the model object now
        model = Model()
        model.setLogger(Logger.IJ_LOGGER)

        # Prepare settings object
        settings = Settings(imp)
        settings.detectorFactory = LogDetectorFactory()
        settings.detectorSettings = {
            'DO_SUBPIXEL_LOCALIZATION': True,
            'RADIUS': float(channel_params[i]['RADIUS']),
            'TARGET_CHANNEL': 1,
            'THRESHOLD': float(channel_params[i]['THRESHOLD']),
            'DO_MEDIAN_FILTERING': False,
        }

        # Configure spot filters - Classical filter on quality
        filter1 = FeatureFilter('QUALITY', float(0.0), True)
        settings.addSpotFilter(filter1)

        # Configure tracker - We want to allow merges and fusions
        settings.trackerFactory = SparseLAPTrackerFactory()
        settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
        settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = False
        settings.trackerSettings['ALLOW_TRACK_MERGING'] = False
        settings.trackerSettings['LINKING_MAX_DISTANCE'] = float(7)
        settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = float(7)
        settings.trackerSettings['MAX_FRAME_GAP'] = int(8)

        # Add ALL the feature analyzers known to TrackMate
        settings.addAllAnalyzers()

        # Configure track filters
        filter2 = FeatureFilter('TRACK_DISPLACEMENT', float(0.0), True)
        settings.addTrackFilter(filter2)

        # Instantiate plugin
        trackmate = TrackMate(model, settings)

        # Process
        ok = trackmate.checkInput()
        if not ok:
            print("Error in checkInput: " + trackmate.getErrorMessage())
            continue

        ok = trackmate.process()
        if not ok:
            print("Error in process: " + trackmate.getErrorMessage())
            continue

        # Display results
        selectionModel = SelectionModel(model)

        ds = DisplaySettingsIO.readUserDefault()
        ds.setTrackColorBy(TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX)
        ds.setSpotColorBy(TrackMateObject.TRACKS, TrackIndexAnalyzer.TRACK_INDEX)

        displayer = HyperStackDisplayer(model, selectionModel, imp, ds)
        displayer.render()
        displayer.refresh()

        # Echo results with the logger we set at start:
        model.getLogger().log(str(model))

        # Export results to XML
        outputXmlPath = channel_params[i]['OUTPUT_XML']
        outFile = File(outputXmlPath)
        ExportTracksToXML.export(model, settings, outFile)
        print("Tracks successfully exported to " + outputXmlPath)

        # Export statistics to CSV
        outputCsvPath = channel_params[i]['OUTPUT_CSV']
        with open(outputCsvPath, 'w') as f:
            f.write("ID,TRACK_ID,QUALITY,POSITION_X,POSITION_Y,POSITION_Z,POSITION_T,FRAME,RADIUS,VISIBILITY,MANUAL_SPOT_COLOR,MEAN_INTENSITY_CH1,MEDIAN_INTENSITY_CH1,MIN_INTENSITY_CH1,MAX_INTENSITY_CH1,TOTAL_INTENSITY_CH1,STD_INTENSITY_CH1,CONTRAST_CH1,SNR_CH1\n")
            track_model = model.getTrackModel()
            if track_model.nTracks(True) == 0:
                print("No spots detected, creating empty output files.")
            else:
                for track_id in track_model.trackIDs(True):
                    track = track_model.trackSpots(track_id)
                    for spot in track:
                        f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
                            spot.ID(),
                            track_id,
                            spot.getFeature('QUALITY'),
                            spot.getDoublePosition(0),
                            spot.getDoublePosition(1),
                            spot.getDoublePosition(2),
                            spot.getFeature('POSITION_T'),
                            spot.getFeature('FRAME'),
                            spot.getFeature('RADIUS'),
                            spot.getFeature('VISIBILITY'),
                            spot.getFeature('MANUAL_SPOT_COLOR'),
                            spot.getFeature('MEAN_INTENSITY_CH1'),
                            spot.getFeature('MEDIAN_INTENSITY_CH1'),
                            spot.getFeature('MIN_INTENSITY_CH1'),
                            spot.getFeature('MAX_INTENSITY_CH1'),
                            spot.getFeature('TOTAL_INTENSITY_CH1'),
                            spot.getFeature('STD_INTENSITY_CH1'),
                            spot.getFeature('CONTRAST_CH1'),
                            spot.getFeature('SNR_CH1')
                        ))
        imp.close()

# Main script
if __name__ == "__main__":
    if input_dir is not None:
        nd_files = [f for f in os.listdir(input_dir) if f.endswith('.nd')]
        for nd_file in nd_files:
            process_nd_file(os.path.join(input_dir, nd_file))
        print("Analysis complete.")
        #IJ.run("Quit")  # Close ImageJ when done
    else:
        print("No directory selected.")
