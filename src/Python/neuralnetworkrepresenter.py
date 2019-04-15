import minieigen
import numpy as np

class NeuralNetworkRepresenter(object):
	""" Class acts as interface object between the neural network's loss given by 
	TATi and and the function to minimize as required by BASSO.
	"""

	def __init__(self, nn, dataset_dimension, true_labels):
		""" Cstor of class.
		
		Args:
		  nn: neural network class
		  dataset_dimension: dimension of dataset
		  true_labels: set of true labels of the dataset
		"""
		self.output_dimension = nn.get_options().output_dimension
		self.numpy_params = np.zeros(nn.num_parameters())
		self.labels = np.zeros(dataset_dimension*self.output_dimension)
		self.grads = np.zeros((dataset_dimension*self.output_dimension,nn.num_parameters()))
		self.nn = nn
		self.dataset_dimension = dataset_dimension
		self.true_labels = true_labels

	def _convertEigenToNumpy(self, _argument):
		""" Convert parameters from Eigen::VectorXd to numpy.ndarray.
		Args:
		   _argument: minieigen array to convert
		   
		Returns:
		   argument converted to numpy
		"""
		for i in range(self.nn.num_parameters()):
			self.numpy_params[i] = _argument[i]
		return self.numpy_params

	def mapping(self, _argument):
		""" Maps given _argument in the feature space to the label space.
		
		Args:
		  argument: features as numpy array 
		  
		Returns:
		  labels predicted by neural network
		"""
		print("Calculating mapping at "+str(_argument))
		# assign parameters
		self.nn.parameters = self._convertEigenToNumpy(_argument)
		# get predicted labels for the whole dataset
		for i in range(self.dataset_dimension):
			self.nn.loss()
			#print(self.nn.dataset)
			temp_label = self.nn.predict(features=np.asarray(self.nn.dataset[0]))[0][0]
			#print(temp_label)
			self.labels[i*self.output_dimension:(i+1)*self.output_dimension] = temp_label
		print("Predicted labels: "+str(self.labels))
		print("Residual: "+str(self.labels-np.reshape(self.true_labels, (self.dataset_dimension*self.output_dimension))))
		return minieigen.VectorX(self.labels)

	def derivative(self, _argument):
		""" Jacobian of the mapping function at _argument.
		
		Args:
		  _argument: features as numpy array 

		Returns:
		  Jacobian at _argument
		"""
		print("Calculating jacobian at "+str(_argument))
		# assign parameters
		self.nn.parameters = self._convertEigenToNumpy(_argument)
		# return gradients
		for i in range(self.dataset_dimension):
			temp_grads = self.nn.gradients()
			#print(temp_grads)
			temp_label = self.nn.predict(features=np.asarray(self.nn.dataset[0]))[0]
			for j in range(self.output_dimension):
				# gradient of representer is hidden inside gradient of loss
				# factor of 2 depends on mean squared loss
				print(i,j,self.grads.shape, self.true_labels.shape,temp_label.shape,temp_grads.shape)
				difference = self.true_labels[i,j] - temp_label[0,j]
				if difference != 0:
					self.grads[i*self.output_dimension+j,:] = -1./(2.*difference)*temp_grads
				else:
					self.grads[i*self.output_dimension+j,:] = 0.*temp_grads
				#print(self.nn.dataset)
		print("Jacobian: "+str(self.grads))
		return self.grads.T

