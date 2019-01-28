cat("mb.c.model(model)", "\n")

mb.c.model <- function(model){

	xb.hat<-model$linear

	n<-length(xb.hat)
	ord<-order(xb.hat)
	xb.hat.ord<-xb.hat[ord]
	p.hat<-plogis(xb.hat.ord)
	q.hat<-1-p.hat
	V1<-(p.hat*(cumsum(q.hat)-q.hat)+q.hat*(sum(p.hat)-cumsum(p.hat)))/(n-1)
	V2<-(p.hat*(sum(q.hat)-q.hat)+q.hat*(sum(p.hat)-p.hat))/(n-1)

	U1<-sum(V1)/n
	U2<-sum(V2)/n
	mb.c<-U1/U2

	s11<-sum((V1-U1)*(V1-U1))/(n-1)
	s22<-sum((V2-U2)*(V2-U2))/(n-1)
	s12<-sum((V1-U1)*(V2-U2))/(n-1)
	s_sq<-4*(U2*U2*s11-2*U1*U2*s12+U1*U1*s22)/(U2^4)
	se<-sqrt(s_sq/n)


	beta<-model$coef
	var.beta<-model$var
	se.beta<-sqrt(diag(var.beta))
	x<-cbind(1,model$x[ord,])
	m<-length(beta)
	D<-rep(0,m)

	for (i in 1:m)
		{
		p.hat<-plogis(xb.hat.ord+se.beta[i]*x[,i])
		q.hat<-1-p.hat
		V1<-(p.hat*(cumsum(q.hat)-q.hat)+q.hat*(sum(p.hat)-cumsum(p.hat)))/(n-1)
		V2<-(p.hat*(sum(q.hat)-q.hat)+q.hat*(sum(p.hat)-p.hat))/(n-1)
		U1<-sum(V1)/n
		U2<-sum(V2)/n
		mb.c.plus<-U1/U2

		p.hat<-plogis(xb.hat.ord-se.beta[i]*x[,i])
		q.hat<-1-p.hat
		V1<-(p.hat*(cumsum(q.hat)-q.hat)+q.hat*(sum(p.hat)-cumsum(p.hat)))/(n-1)
		V2<-(p.hat*(sum(q.hat)-q.hat)+q.hat*(sum(p.hat)-p.hat))/(n-1)
		U1<-sum(V1)/n
		U2<-sum(V2)/n
		mb.c.min<-U1/U2
	
		D[i]<-(mb.c.plus-mb.c.min)/(2*se.beta[i])
		}

	DT.varbeta.D<-t(D)%*%var.beta%*%D

	return(list(mbc=mb.c,se=sqrt(s_sq/n+DT.varbeta.D)))
}

