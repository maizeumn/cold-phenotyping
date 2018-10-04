import argparse
import os
import sys
import cv2
import numpy as np
from matplotlib import pyplot as plt
from plantcv import plantcv as pcv
from plantcv.plantcv import params

def options():
	parser = argparse.ArgumentParser(description="Imaging processing with PlantCV.")
	parser.add_argument("-d", "--image", help="Input image file.", required=True)
	parser.add_argument("-r","--result", help="Result file.", required=False)
	parser.add_argument("-p", "--pdfs", help="Naive Bayes PDF file.", required=True)
	parser.add_argument("-o", "--outdir", help="Output directory for image files.", required=False)
	parser.add_argument("-w","--writeimg", help="Write out images.", default=False, action="store_true")
	parser.add_argument("-D", "--debug", help="Turn on debug, prints intermediate images.", action=None)
	args = parser.parse_args()
	return args

def main():
	args = options()
	
	os.chdir(args.outdir)
	
	# Read RGB image
	img, path, filename = pcv.readimage(args.image, mode="native")

	# Get metadata from file name
	geno_name = filename.split("}{")
	geno_name = geno_name[5]
	geno_name = geno_name.split("_")
	geno_name = geno_name[1]

	day = filename.split("}{")
	day = day[7]
	day = day.split("_")
	day= day[1]
	day = day.split("}")
	day = day[0]

	plot = filename.split("}{")
	plot = plot[0]
	plot = plot.split("_")
	plot = plot[1]

	exp_name = filename.split("}{")
	exp_name = exp_name[1]
	exp_name = exp_name.split("_")
	exp_name = exp_name[1]

	treat_name = filename.split("}{")
	treat_name = treat_name[6]

	# Create masks using Naive Bayes Classifier and PDFs file
	masks = pcv.naive_bayes_classifier(img, args.pdfs)

	# The following code will identify the racks in the image, find the top edge, and choose a line along the edge to pick a y coordinate to trim any soil/pot pixels identified as plant material.
	
	# Convert RGB to HSV and extract the Value channel
	v = pcv.rgb2gray_hsv(img, 'v')
	
	# Threshold the Value image
	v_thresh = pcv.threshold.binary(v, 98, 255, 'light')
	
	# Dilate mask to fill holes
	dilate_racks = pcv.dilate(v_thresh, 2, 1)
	
	# Fill in small objects
	mask = np.copy(dilate_racks)
	fill_racks = pcv.fill(mask, 100000)

	#edge detection
	edges = cv2.Canny(fill_racks,60,180)
	
	#write all the straight lines from edge detection
	lines = cv2.HoughLinesP(edges,rho = 1,theta = 1*np.pi/180,threshold = 150,minLineLength = 50,maxLineGap = 15)
	N = lines.shape[0]
	for i in range(N):
		x1 = lines[i][0][0]
		y1 = lines[i][0][1]
		x2 = lines[i][0][2]
		y2 = lines[i][0][3]
		cv2.line(img,(x1,y1),(x2,y2),(255,0,0),2)

	# keep only horizontal lines
	N = lines.shape[0]
	tokeep = []
	for i in range(N):
		want = (abs(lines[i][0][1]-lines[i][0][3]))<=10
		tokeep.append(want)

	lines = lines[tokeep]

	# keep only lines in lower half of image
	N = lines.shape[0]
	tokeep = []
	for i in range(N):
		want = 3100>lines[i][0][1]>2300
		tokeep.append(want)

	lines = lines[tokeep]

	# assign lines to positions around plants
	N = lines.shape[0]
	tokeep = []
	left = []
	mid = []
	right = []

	for i in range(N):
		leftones = lines[i][0][2]<=2000
		left.append(leftones)
		
		midones = 3000>lines[i][0][2]>2000
		mid.append(midones)
		
		rightones = lines[i][0][0]>=3300
		right.append(rightones)
		
	right = lines[right]
	left = lines[left]
	mid = lines[mid]
	
	# choose y values for right left mid adding some pixels to go about the pot (subtract because of orientation of axis)
	y_left = left[0][0][3]-50
	y_mid = mid[0][0][3]-50
	y_right= right[0][0][3]-50
	
	# reload original image to write new lines on
	img, path, filename = pcv.readimage(args.image)

	# write horizontal lines on image
	cv2.line(img,(left[0][0][0],left[0][0][1]),(left[0][0][2],left[0][0][3]),(255,255,51),2)
	cv2.line(img,(mid[0][0][0],mid[0][0][1]),(mid[0][0][2],mid[0][0][3]),(255,255,51),2)  
	cv2.line(img,(right[0][0][0],right[0][0][1]),(right[0][0][2],right[0][0][3]),(255,255,51),2)
	
	# Add masks together
	added = masks["healthy"] + masks["necrosis"] + masks["stem"]
	
	# Dilate mask to fill holes
	dilate_img = pcv.dilate(added, 2, 1)
	
	# Fill in small objects
	mask = np.copy(dilate_img)
	fill_img = pcv.fill(mask, 400)
	
	ret, inverted = cv2.threshold(fill_img, 75, 255, cv2.THRESH_BINARY_INV)
	
	# Dilate mask to fill holes of plant
	dilate_inv = pcv.dilate(inverted, 2, 1)
	
	# Fill in small objects of plant
	mask2 = np.copy(dilate_inv)
	fill_plant = pcv.fill(mask2, 20)
	
	inverted_img = pcv.invert(fill_plant)
	
	# Identify objects
	id_objects, obj_hierarchy = pcv.find_objects(img, inverted_img)
	
	# Define ROIs
	roi_left, roi_hierarchy_left = pcv.roi.rectangle(280, 1280, 1275, 1200, img)
	roi_mid, roi_hierarchy_mid = pcv.roi.rectangle(1900, 1280, 1275, 1200, img)
	roi_right, roi_hierarchy_right = pcv.roi.rectangle(3600, 1280, 1275, 1200, img)
	
	# Decide which objects to keep
	roi_objects_left, roi_obj_hierarchy_left, kept_mask_left, obj_area_left = pcv.roi_objects(img, 'partial', roi_left, roi_hierarchy_left,id_objects, obj_hierarchy)
	roi_objects_mid, roi_obj_hierarchy_mid, kept_mask_mid, obj_area_mid = pcv.roi_objects(img, 'partial', roi_mid, roi_hierarchy_mid,id_objects, obj_hierarchy)
	roi_objects_right, roi_obj_hierarchy_right, kept_mask_right, obj_area_right = pcv.roi_objects(img, 'partial', roi_right, roi_hierarchy_right,id_objects, obj_hierarchy)
	
	# Combine objects
	obj_r, mask_r = pcv.object_composition(img, roi_objects_right, roi_obj_hierarchy_right)
	obj_m, mask_m = pcv.object_composition(img, roi_objects_mid, roi_obj_hierarchy_mid)
	obj_l, mask_l = pcv.object_composition(img, roi_objects_left, roi_obj_hierarchy_left)
	
	def analyze_bound_horizontal2(img, obj, mask, line_position, filename=False):
		
		ori_img = np.copy(img)

		# Draw line horizontal line through bottom of image, that is adjusted to user input height
		if len(np.shape(ori_img)) == 3:
			iy, ix, iz = np.shape(ori_img)
		else:
			iy, ix = np.shape(ori_img)
		size = (iy, ix)
		size1 = (iy, ix, 3)
		background = np.zeros(size, dtype=np.uint8)
		wback = np.zeros(size1, dtype=np.uint8)
		x_coor = int(ix)
		y_coor = int(iy) - int(line_position)
		rec_corner = int(iy - 2)
		rec_point1 = (1, rec_corner)
		rec_point2 = (x_coor - 2, y_coor - 2)
		cv2.rectangle(background, rec_point1, rec_point2, (255), 1)
		below_contour, below_hierarchy = cv2.findContours(background, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)[-2:]

		below = []
		above = []
		mask_nonzerox, mask_nonzeroy = np.nonzero(mask)
		obj_points = np.vstack((mask_nonzeroy, mask_nonzerox))
		obj_points1 = np.transpose(obj_points)

		for i, c in enumerate(obj_points1):
			xy = tuple(c)
			pptest = cv2.pointPolygonTest(below_contour[0], xy, measureDist=False)
			if pptest == 1:
				below.append(xy)
				cv2.circle(ori_img, xy, 1, (0, 0, 255))
				cv2.circle(wback, xy, 1, (0, 0, 0))
			else:
				above.append(xy)
				cv2.circle(ori_img, xy, 1, (0, 255, 0))
				cv2.circle(wback, xy, 1, (255, 255, 255))

		return wback
	
	ori_img = np.copy(img)

	# Draw line horizontal line through bottom of image, that is adjusted to user input height
	if len(np.shape(ori_img)) == 3:
		iy, ix, iz = np.shape(ori_img)
	else:
		iy, ix = np.shape(ori_img)
	
	if obj_r is not None:
		wback_r = analyze_bound_horizontal2(img, obj_r, mask_r, iy-y_right)
	if obj_m is not None:
		wback_m = analyze_bound_horizontal2(img, obj_m, mask_m, iy-y_mid)
	if obj_l is not None:
		wback_l = analyze_bound_horizontal2(img, obj_l, mask_l, iy-y_left)

	threshold_light = pcv.threshold.binary(img, 1, 1, 'dark')

	if obj_r is not None:
		fgmask_r = pcv.background_subtraction(wback_r, threshold_light)
	if obj_m is not None:
		fgmask_m = pcv.background_subtraction(wback_m, threshold_light)
	if obj_l is not None:
		fgmask_l = pcv.background_subtraction(wback_l, threshold_light)

	if obj_l is not None:
		id_objects_left, obj_hierarchy_left = pcv.find_objects(img, fgmask_l)
	if obj_m is not None:
		id_objects_mid, obj_hierarchy_mid = pcv.find_objects(img, fgmask_m)
	if obj_r is not None:
		id_objects_right, obj_hierarchy_right = pcv.find_objects(img, fgmask_r)

	# Combine objects
	if obj_r is not None:
		obj_r2, mask_r2 = pcv.object_composition(img, id_objects_right, obj_hierarchy_right)
	if obj_m is not None:
		obj_m2, mask_m2 = pcv.object_composition(img, id_objects_mid, obj_hierarchy_mid)
	if obj_l is not None:
		obj_l2, mask_l2 = pcv.object_composition(img, id_objects_left, obj_hierarchy_left)

	# Shape measurements
	if obj_l is not None:
		shape_header_left, shape_data_left, shape_img_left = pcv.analyze_object(img, obj_l2, fgmask_l, geno_name+'_'+plot+'_'+'A'+'_'+day+'_'+'shape.jpg')

	if obj_r is not None:
		shape_header_right, shape_data_right, shape_img_right = pcv.analyze_object(img, obj_r2, fgmask_r, geno_name+'_'+plot+'_'+'C'+'_'+day+'_'+'shape.jpg')

	if obj_m is not None:
		shape_header_mid, shape_data_mid, shape_img_mid = pcv.analyze_object(img, obj_m2, fgmask_m, geno_name+'_'+plot+'_'+'B'+'_'+day+'_'+'shape.jpg')

	# Color data
	if obj_r is not None:
		color_header_right, color_data_right, norm_slice_right = pcv.analyze_color(img, fgmask_r, 256, None,'v','img', geno_name+'_'+plot+'_'+'C'+'_'+day+'_'+'color.jpg')

	if obj_m is not None:
		color_header_mid, color_data_mid, norm_slice_mid = pcv.analyze_color(img, fgmask_m, 256, None,'v','img', geno_name+'_'+plot+'_'+'B'+'_'+day+'_'+'color.jpg')

	if obj_l is not None:
		color_header_left, color_data_left, norm_slice_left = pcv.analyze_color(img, fgmask_l, 256, None,'v','img', geno_name+'_'+plot+'_'+'A'+'_'+day+'_'+'color.jpg')

	new_header = ['experiment', 'day', 'genotype', 'treatment', 'plot', 'plant', 'percent.necrosis', 'area', 'hull-area', 'solidity', 'perimeter', 'width', 'height', 'longest_axis', 'center-of-mass-x', 'center-of-mass-y', 'hull_vertices', 'in_bounds', 'ellipse_center_x', 'ellipse_center_y', 'ellipse_major_axis', 'ellipse_minor_axis', 'ellipse_angle', 'ellipse_eccentricity', 'bin-number', 'bin-values', 'blue', 'green', 'red', 'lightness', 'green-magenta', 'blue-yellow', 'hue', 'saturation', 'value']
	table = []
	table.append(new_header)

	added2 = masks["healthy"] + masks["stem"]
	
	# Object combine kept objects
	if obj_l is not None:
		masked_image_healthy_left = pcv.apply_mask(added2, fgmask_l, 'black')
		masked_image_necrosis_left = pcv.apply_mask(masks["necrosis"], fgmask_l, 'black')
		added_obj_left = masked_image_healthy_left + masked_image_necrosis_left

		sample = "A"
		
		# Calculations
		necrosis_left = np.sum(masked_image_necrosis_left)
		necrosis_percent_left = float(necrosis_left)/np.sum(added_obj_left)
		
		table.append([exp_name, day, geno_name, treat_name, plot, sample, round(necrosis_percent_left, 5), shape_data_left[1], shape_data_left[2], shape_data_left[3], shape_data_left[4], shape_data_left[5], shape_data_left[6], shape_data_left[7], shape_data_left[8], shape_data_left[9], shape_data_left[10], shape_data_left[11], shape_data_left[12], shape_data_left[13], shape_data_left[14], shape_data_left[15], shape_data_left[16], shape_data_left[17], '"{}"'.format(color_data_left[1]), '"{}"'.format(color_data_left[2]), '"{}"'.format(color_data_left[3]), '"{}"'.format(color_data_left[4]), '"{}"'.format(color_data_left[5]), '"{}"'.format(color_data_left[6]), '"{}"'.format(color_data_left[7]), '"{}"'.format(color_data_left[8]), '"{}"'.format(color_data_left[9]), '"{}"'.format(color_data_left[10]), '"{}"'.format(color_data_left[11])])

	# Object combine kept objects
	if obj_m is not None:
		masked_image_healthy_mid = pcv.apply_mask(added2, fgmask_m, 'black')
		masked_image_necrosis_mid = pcv.apply_mask(masks["necrosis"], fgmask_m, 'black')
		added_obj_mid = masked_image_healthy_mid + masked_image_necrosis_mid
		
		sample = "B"
		
		# Calculations
		necrosis_mid = np.sum(masked_image_necrosis_mid)
		necrosis_percent_mid = float(necrosis_mid)/np.sum(added_obj_mid)
		
		table.append([exp_name, day, geno_name, treat_name, plot, sample,  round(necrosis_percent_mid, 5), shape_data_mid[1], shape_data_mid[2], shape_data_mid[3], shape_data_mid[4], shape_data_mid[5], shape_data_mid[6], shape_data_mid[7], shape_data_mid[8], shape_data_mid[9], shape_data_mid[10], shape_data_mid[11], shape_data_mid[12], shape_data_mid[13], shape_data_mid[14], shape_data_mid[15], shape_data_mid[16], shape_data_mid[17], '"{}"'.format(color_data_mid[1]), '"{}"'.format(color_data_mid[2]), '"{}"'.format(color_data_mid[3]), '"{}"'.format(color_data_mid[4]), '"{}"'.format(color_data_mid[5]), '"{}"'.format(color_data_mid[6]), '"{}"'.format(color_data_mid[7]), '"{}"'.format(color_data_mid[8]), '"{}"'.format(color_data_mid[9]), '"{}"'.format(color_data_mid[10]), '"{}"'.format(color_data_mid[11])])

	# Object combine kept objects
	if obj_r is not None:
		masked_image_healthy_right = pcv.apply_mask(added2, fgmask_r, 'black')
		masked_image_necrosis_right = pcv.apply_mask(masks["necrosis"], fgmask_r, 'black')
		added_obj_right = masked_image_healthy_right + masked_image_necrosis_right
		
		sample = "C"
		
		# Calculations
		necrosis_right = np.sum(masked_image_necrosis_right)
		necrosis_percent_right = float(necrosis_right)/np.sum(added_obj_right)
		
		table.append([exp_name, day, geno_name, treat_name, plot, sample,  round(necrosis_percent_right, 5), shape_data_right[1], shape_data_right[2], shape_data_right[3], shape_data_right[4], shape_data_right[5], shape_data_right[6], shape_data_right[7], shape_data_right[8], shape_data_right[9], shape_data_right[10], shape_data_right[11], shape_data_right[12], shape_data_right[13], shape_data_right[14], shape_data_right[15], shape_data_right[16], shape_data_right[17], '"{}"'.format(color_data_right[1]), '"{}"'.format(color_data_right[2]), '"{}"'.format(color_data_right[3]), '"{}"'.format(color_data_right[4]), '"{}"'.format(color_data_right[5]), '"{}"'.format(color_data_right[6]), '"{}"'.format(color_data_right[7]), '"{}"'.format(color_data_right[8]), '"{}"'.format(color_data_right[9]), '"{}"'.format(color_data_right[10]), '"{}"'.format(color_data_right[11])])

	if obj_l is not None:	
		merged2 = cv2.merge([masked_image_healthy_left, np.zeros(np.shape(masks["healthy"]),dtype=np.uint8), masked_image_necrosis_left]) #blue, green, red
		pcv.print_image(merged2, geno_name+'_'+plot+'_'+'A'+'_'+day+'_'+'merged.jpg')
	if obj_m is not None:
		merged3 = cv2.merge([masked_image_healthy_mid, np.zeros(np.shape(masks["healthy"]),dtype=np.uint8), masked_image_necrosis_mid]) #blue, green, red
		pcv.print_image(merged3, geno_name+'_'+plot+'_'+'B'+'_'+day+'_'+'merged.jpg')
	if obj_r is not None:
		merged4 = cv2.merge([masked_image_healthy_right, np.zeros(np.shape(masks["healthy"]),dtype=np.uint8), masked_image_necrosis_right]) #blue, green, red
		pcv.print_image(merged4, geno_name+'_'+plot+'_'+'C'+'_'+day+'_'+'merged.jpg')
	
	# Save area results to file (individual csv files for one image...)
	file_name = filename.split("}{")
	file_name = file_name[0]+"}{"+file_name[5]+"}{"+file_name[7]

	outfile = str(file_name[:-4]) + 'csv' 
	with open(outfile, 'w') as f:
		for row in table:
			f.write(','.join(map(str, row)) + '\n')
			
	print(filename)

if __name__ == '__main__':
	main()